#include <mpi.h>
#include <stdio.h>
int main (int argc, char *argv [])
{
	for (int i = 0; i < argc; i++)
	{
		printf ("%s\n", argv [i]);
	}

	int commsize, my_rank;
	printf ("argc before calling MPI_Init() = %d\n", argc);
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &commsize);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	printf ("Communicator size = %d My rank= %d\n", commsize, my_rank);
	printf ("argc after calling MPI_Init() = %d\n", argc);
	MPI_Finalize ();
}