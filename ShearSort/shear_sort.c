#define _CRT_SECURE_NO_WARNINGS
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MASTER 0
#define MAX_WORD_LEN 20
#define NUM_OF_DIMENSIONS 2
#define ROW 0
#define COL 1
#define ASCENDING 0
#define DESCENDING 1


//prototypesmentation
char* read_words_from_file(const char* file_path, int* num_of_words);
void print_words(char* words, int num_of_words);
void check_allocation(void* ptr);
void set_matrix_dimensions(int* dimensions, int num_of_words);
void validate_num_of_procs(int num_of_procs, int num_of_words);
int ascending(char* word, char* other_word);
int descending(char* word, char* other_word);
void print_matrix(char* words, int dimension, int size);
double log_base2(double num);
void odd_even_sort(char* value, int(*compare)(char*, char*), int firstRank, int secondRank, int location, int size, MPI_Comm comm);
void shear_sort(char* value, MPI_Comm cart_comm, int size);





int main(int argc, char *argv[])
{
	//change file path to the text file
	const char* file_path = "C:\\Users\\Ran\\Desktop\\InitialMPIProject\\words.txt";

	char* words = NULL, *result = NULL;
	char word[MAX_WORD_LEN];
	int dimensions[NUM_OF_DIMENSIONS], period[NUM_OF_DIMENSIONS] = { 0 }, reorder = 0;
	int  numprocs, rank, num_of_words;
	MPI_Comm cart_comm;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);


	//master read from txt file and print it
	if (rank == MASTER)
	{
		words = read_words_from_file(file_path, &num_of_words);
		validate_num_of_procs(numprocs, num_of_words);

		result = (char*)malloc(sizeof(char)* MAX_WORD_LEN * (num_of_words));
		check_allocation(result);
	}

	//send the number of words to everyone (n*n size matrix)
	MPI_Bcast(&num_of_words, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

	//set dimensions
	set_matrix_dimensions(dimensions, num_of_words);

	//master send each process one word
	MPI_Scatter(words, MAX_WORD_LEN, MPI_CHAR, word, MAX_WORD_LEN, MPI_CHAR, MASTER, MPI_COMM_WORLD);

	//create cartesian topology 
	MPI_Cart_create(MPI_COMM_WORLD, NUM_OF_DIMENSIONS, dimensions, period, reorder, &cart_comm);

	//sort by lexicographic order
	shear_sort(word, cart_comm, dimensions[0]);

	//master recieve data from each process
	MPI_Gather(word, MAX_WORD_LEN, MPI_CHAR, result, MAX_WORD_LEN, MPI_CHAR, MASTER, cart_comm);

	//master print the result - snake-like sorted matrix
	if (rank == MASTER)
	{
		printf("\n\n////////////////////////  UNSORTED  ///////////////////////\n");
		printf("//////////////////////////////////////////////////////////\n\n");
		print_matrix(words, dimensions[0], num_of_words);
		printf("\n\n////////////////////////  SORTED  /////////////////////////\n");
		printf("//////////////////////////////////////////////////////////\n");
		print_matrix(result, dimensions[0], num_of_words);
	}

	MPI_Finalize();
	return 0;
}


char* read_words_from_file(const char* file_path, int* num_of_words)
{
	FILE* file;
	char* words;
	int num, i;

	//open file to read
	file = fopen(file_path, "r");

	//file not opened
	if (file == NULL)
	{
		printf("file not opened - make sure the .txt file is in the correct path");
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	//read first line 
	fscanf(file, "%d", &num);

	//set the total number of words
	*num_of_words = num;

	//allocate and read values from text to linear array 
	words = (char*)malloc(sizeof(char)* MAX_WORD_LEN * (*num_of_words));
	check_allocation(words);



	for (i = 0; i < *num_of_words; i++)
	{
		fscanf(file, "%s", words + i*MAX_WORD_LEN);
	}
	return words;
}


void print_words(char* words, int num_of_words)
{
	int i;
	char* ptr = words;

	for (i = 0; i < num_of_words; i++)
	{
		printf("%.*s\n", MAX_WORD_LEN, ptr);
		ptr += MAX_WORD_LEN;
	}
}


void set_matrix_dimensions(int* dimensions, int num_of_words)
{
	double square;

	//check if the num of words in the text file is valid number
	square = sqrt(num_of_words);

	if (square != rint(square))
	{
		printf("incorrect num of words");
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	//set the column & row to the square number
	dimensions[0] = dimensions[1] = (int)square;
}


void check_allocation(void* ptr)
{
	if (!ptr)
	{
		printf("No memory for allocation");
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}


void validate_num_of_procs(int num_of_procs, int num_of_words)
{
	//necessary number of procs is not equal to the actual number
	if (num_of_procs != num_of_words)
	{
		printf("number of proccess need to be as the number of words = %d", num_of_words);
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

//shear sort algorithm using odd even sort - assume that n is even
void shear_sort(char* value, MPI_Comm cart_comm, int size)
{
	int num_of_iterations, rank;
	int source, dest, coords[NUM_OF_DIMENSIONS];

	//number of iteration in parallel shear sprt 2logn+1
	num_of_iterations = 2 * ((int)log2(size)) + 1;

	MPI_Comm_rank(cart_comm, &rank);
	MPI_Cart_coords(cart_comm, rank, NUM_OF_DIMENSIONS, coords);

	for (int i = 0; i < num_of_iterations; i++)
	{
		//even phase - sorting rows
		if (i % 2 == 0)
		{
			//get correct neighboors for communication
			MPI_Cart_shift(cart_comm, 1, 1, &source, &dest);
			if (coords[ROW] % 2 == 0)
			{
				odd_even_sort(value, ascending, source, dest, coords[COL], size, cart_comm);
			}
			else
			{
				odd_even_sort(value, descending, source, dest, coords[COL], size, cart_comm);
			}
		}
		//odd phase - sorting columns
		else
		{
			//get correct neighboors for communication
			MPI_Cart_shift(cart_comm, 0, 1, &source, &dest);
			odd_even_sort(value, ascending, source, dest, coords[ROW], size, cart_comm);
		}
	}
}


void odd_even_sort(char* value, int(*compare)(char*, char*), int source, int dest, int location, int size, MPI_Comm comm) {
	char neighboor_value[MAX_WORD_LEN];
	MPI_Status status;

	for (int i = 0; i < size; i++)
	{
		//sending
		if ((i % 2 == location % 2) && (location < size - 1))
		{
			//send value to neighboor and update by the correct value from comare function
			MPI_Send(value, MAX_WORD_LEN, MPI_CHAR, dest, 0, comm);
			MPI_Recv(value, MAX_WORD_LEN, MPI_CHAR, dest, 0, comm, &status);
		}
		else if ((i % 2 != location % 2) && (location > 0))
		{
			//recieving value from neeighboor and decide what to do with it
			MPI_Recv(neighboor_value, MAX_WORD_LEN, MPI_CHAR, source, 0, comm, &status);

			//compare two values according the orientation
			if (compare(neighboor_value, value))
			{
				//swap needed - send my value to neighboor and update mine
				MPI_Send(value, MAX_WORD_LEN, MPI_CHAR, source, 0, comm);
				strcpy(value, neighboor_value);
			}
			else
			{
				//no need for swap - send neighboor his value 
				MPI_Send(neighboor_value, MAX_WORD_LEN, MPI_CHAR, source, 0, comm);
			}
		}
	}
}


int ascending(char* word, char* other_word)
{
	return strcmp(word, other_word) > 0;
}


int descending(char* word, char* other_word)
{
	return strcmp(other_word, word) > 0;
}


void print_matrix(char* words, int dimension, int size)
{
	int i;
	char* ptr = words;

	for (i = 0; i < size; i++)
	{
		if (i % dimension == 0)
			printf("\n");

		printf("%-15.*s", MAX_WORD_LEN, ptr);
		ptr += MAX_WORD_LEN;
	}
	printf("\n\n");
}


double log_base2(double num)
{
	return ((log(num)) / (log(2.0)));
}