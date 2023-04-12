#include <mpi.h>
#include <stdio.h>

int main (int argc, char *argv [])
{
	int my_rank;
	MPI_Status status;
	unsigned int num = 0;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

	if (my_rank == 0)
	{
		double t_start = MPI_Wtime ();

		num = 10;
		MPI_Send (&num, 1, MPI_UNSIGNED, 1, 0, MPI_COMM_WORLD);
		
		double t_stop = MPI_Wtime ();
		printf ("time = %f s\n", t_stop - t_start);
	}
	else if (my_rank == 1)
	{
		MPI_Recv (&num, 1, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		printf ("Received num = %d\n", num);
	}

	MPI_Finalize ();
	return 0;
}