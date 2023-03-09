#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

double sum_elements (int my_rank, unsigned int epp, int N);

int main (int argc, char *argv [])
{
	MPI_Init (&argc, &argv);

	int commsize, my_rank;
	MPI_Status status;
	MPI_Comm_size (MPI_COMM_WORLD, &commsize);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

	int var = 0;

	if (my_rank == 0)
	{
		var = atoi (argv [1]);
		printf ("var = %d (initial value)\n", var);

		var = var++;
		printf ("My rank = %d, var = %d\n", my_rank, var);

		MPI_Send (&var, 1, MPI_INT, my_rank + 1, 0, MPI_COMM_WORLD);
		MPI_Recv (&var, 1, MPI_INT, commsize - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		printf ("My rank = %d, var = %d (received from process #%d)\n", my_rank, var, commsize - 1);
	}
	else
	{
		MPI_Recv (&var, 1, MPI_INT, my_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		var = var++;
		printf ("My rank = %d, var = %d\n", my_rank, var);

		if (my_rank < commsize - 1)
		{
			MPI_Send (&var, 1, MPI_INT, my_rank + 1, 0, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Send (&var, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}

	MPI_Finalize ();
	return 0;
}