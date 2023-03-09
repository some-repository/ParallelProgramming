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

	int N = atoi (argv [1]);
	unsigned int epp = 0;

	if (my_rank == 0)
	{
		//printf ("N = %d\n", N);

		if (N % commsize == 0)
		{
			epp = N / commsize; // elements per process
		}
		else
		{
			epp = (N / commsize) + 1; // elements per process
		}

		for (int i = 1; i < commsize; i++)
		{
			MPI_Send (&epp, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD); // tag = 0
		}

		double sum = sum_elements (my_rank, epp, N);

		double tmp = 0;
		for (int i = 1; i < commsize; i++)
		{
			MPI_Recv (&tmp, 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			sum += tmp;
		}

		printf ("sum = %f\n", sum);
	}
	else
	{
		MPI_Recv (&epp, 1, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		double sum = sum_elements (my_rank, epp, N);

		MPI_Send (&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // tag = 0
	}

	MPI_Finalize ();
	return 0;
}

double sum_elements (int my_rank, unsigned int epp, int N)
{
	double sum = 0;
	for (int i = my_rank * epp + 1; i < (my_rank + 1) * epp + 1; i++)
	{
		if (i > N)
		{
			break;
		}
		sum += 1 / (double) i;
		//printf ("My rank = %d, i = %d\n", my_rank, i);
	}
	return sum;
}