#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/unistd.h>
#include <math.h>
#include <mpi.h>

// solve differential equation U't + a*U'x = f (x, t)

double psi (double t); // U (0, t) = psi (t)
double phi (double x); // U (x, 0) = phi (t)
double f (double x, double t); // the right part of the differential equation
double U_analytical (double x, double t);
double find_max_error (double ** U_arr, size_t str_num, size_t col_num, double tau, double h);

double ** allocate_2D_array (size_t str_num, size_t col_num);
void free_2D_array (double ** array, size_t str_num);

int main (void)
{
	int commsize, my_rank;
	MPI_Status status;
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size (MPI_COMM_WORLD, &commsize);

	const double a = 2; // coefficient of the differential equation

	const double T = 1; // maximum value of t axis
	const int K = 20; // maximum number of element of t axis
	const double tau = T / K; // time step

	const double X = 1; // maximum value of x axis
	const int M = 20; // maximum number of element of x axis
	const double h = X / M; // x step

	// main code should be placed here
	size_t k = my_rank; // number of t moment

	if (my_rank == 0)
	{
		double ** U_arr = allocate_2D_array (K + 1, M + 1);
		for (int m = 0; m < M + 1; m++) // set boundary condition for t = 0
		{
			double x = m * h;
			U_arr [0][m] = phi (x);
		}

		printf ("Max error = %E\n", find_max_error (U_arr, K + 1, M + 1, tau, h));
		print_array_to_file ("solution.csv", U_arr, K + 1, M + 1);		
		free_2D_array (U_arr, K + 1);
	}
	else
	{
		double * x_arr = (double *) calloc (M + 1, sizeof (double *)); // allocate an array for x values in a fixed t moment
    	if (x_arr == NULL)
    	{
    		fprintf (stderr, "x_arr pointer is NULL\n");
    		exit (EXIT_FAILURE);
    	}

    	x_arr [0] = psi (k * tau); // set boundary condition for x = 0

    	free (x_arr);
	}

	MPI_Finalize ();
	return 0;
}

double psi (double t)
{
	return exp (-t);
}

double phi (double x)
{
	return cos (M_PI * x);
}

double f (double x, double t)
{
	return x + t;
}

double U_analytical (double x, double t)
{
	if ((2 * t - x) <= 0)
	{
		return (x * t - (t * t) / 2 + cos (M_PI * (2 * t - x)));
	}
	else
	{
		return (x * t - (t * t) / 2 + (2 * t - x) * (2 * t - x) / 8 + exp (-(t - x / 2)));
	}
}

double find_max_error (double ** U_arr, size_t str_num, size_t col_num, double tau, double h)
{
	double max_error = 0;
	double error = 0;
	for (int k = 0; k < str_num; k++)
	{
		for (int m = 0; m < col_num; m++)
		{
			error = fabs (U_arr [k][m] - U_analytical (m * h, k * tau));
			if (error > max_error)
			{
				max_error = error;
			}
		}
	}

	return max_error;
}

double ** allocate_2D_array (size_t str_num, size_t col_num)
{
	double ** arr = (double **) calloc (str_num, sizeof (double *));
    if (arr == NULL)
    {
    	fprintf (stderr, "U_arr pointer is NULL\n");
    	exit (EXIT_FAILURE);
    }
    for (int k = 0; k < str_num; k++)
    {
        arr [k] = calloc (col_num, sizeof (double));
        if (arr [k] == NULL)
    	{
    		fprintf (stderr, "U_arr pointer is NULL\n");
    		exit (EXIT_FAILURE);
    	}
    }

    return arr;
}

void free_2D_array (double ** array, size_t str_num)
{
	for (int k = 0; k < str_num; k++)
    {
        free (array [k]);
    }
    free (array);
}