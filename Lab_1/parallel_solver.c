//#define PRINT_SOLUTION_TO_FILE // uncomment this line to get solution in .csv file

#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/unistd.h>
#include <math.h>
#include <mpi.h>

#ifndef M_PI // in some versions math.h doesn't contain M_PI constant
    #define M_PI 3.14159265358979323846
#endif

// solve differential equation U't + a*U'x = f (x, t)

double psi (double t); // U (0, t) = psi (t)
double phi (double x); // U (x, 0) = phi (t)
double f (double x, double t); // the right part of the differential equation
double U_analytical (double x, double t);
double find_max_error (double ** U_arr, size_t str_num, size_t col_num, double tau, double h);

double ** allocate_2D_array (size_t str_num, size_t col_num);
void free_2D_array (double ** array, size_t str_num);
void print_array_to_file (const char * filename, double ** arr, size_t str_num, size_t col_num);

int main (int argc, char *argv [])
{
	int commsize, my_rank;
	MPI_Status status;
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size (MPI_COMM_WORLD, &commsize);

	const double a = 2; // coefficient of the differential equation

	const double T = 1; // maximum value of t axis
	const int K = 1000; // maximum number of element of t axis
	const double tau = T / K; // time step

	const double X = 1; // maximum value of x axis
	const int M = 1000; // maximum number of element of x axis
	const double h = X / M; // x step

	size_t task_size = 0;
	size_t common_task_size = 0;
	if ((M + 1) >= commsize)
	{
		common_task_size = (M + 1) / commsize;
	}
	else
	{
		common_task_size = 1;
	}

	if (my_rank == 0)
	{
		double t_start = MPI_Wtime ();

		task_size = common_task_size;

		double ** U_arr = allocate_2D_array (K + 1, M + 1); // here the whole solution will be stored

		for (size_t m = 0; m < task_size; m++) // set boundary condition for t = 0
		{
			double x = m * h;
			U_arr [0][m] = phi (x);
		}
		for (size_t k = 0; k < K + 1; k++) // set boundary condition for x = 0
		{
			double t = k * tau;
			U_arr [k][0] = psi (t);
		}

		for (size_t k = 0; k < K; k++)
		{
			for (size_t m = 1; m < task_size; m++)
			{
				U_arr [k + 1][m] = (2 * f ((m - 0.5) * h, (k + 0.5) * tau) + U_arr [k][m - 1] * ((1 / tau) + (a / h)) + U_arr [k][m] * ((1 / tau) - (a / h)) + U_arr [k + 1][m - 1] * ((a / h) - (1 / tau))) / ((1 / tau) + (a / h));
			}

			MPI_Send (&U_arr [k + 1][task_size - 1], 1, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD);
		}

		size_t proc_used = commsize;
		if (commsize > (M + 1))
		{
			proc_used = M + 1;
		}

		for (size_t m = 1; m < proc_used; m++) // receive parts of solution to the common array U_arr
		{
			size_t count = common_task_size;
			if ((M + 1 - m * common_task_size) < (2 * common_task_size))
			{
				count = common_task_size + ((M + 1) % commsize);
			}

			for (size_t k = 0; k < K + 1; k++)
			{
				MPI_Recv (&U_arr [k][m * common_task_size], count, MPI_DOUBLE_PRECISION, m, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
		}

		double t_stop = MPI_Wtime ();
		printf ("time = %f s\n", t_stop - t_start);
		printf ("Max error = %E\n", find_max_error (U_arr, K + 1, M + 1, tau, h));

		#if defined (PRINT_SOLUTION_TO_FILE)
        	print_array_to_file ("parallel_solution.csv", U_arr, K + 1, M + 1);
    	#endif //PRINT_SOLUTION_TO_FILE
				
		free_2D_array (U_arr, K + 1);
	}
	else
	{
		if ((M + 1) >= commsize)
		{
			if (my_rank == (commsize - 1))
			{
				task_size = common_task_size + ((M + 1) % commsize);
			}
			else
			{
				task_size = common_task_size;
			}
		}
		else
		{
			if (my_rank < (M + 1))
			{
				task_size = 1;
			}
			else
			{
				task_size = 0;
			}
		}

		if (task_size >= 1)
		{
			double ** U_arr_part = allocate_2D_array (K + 1, task_size + 1);

    		for (size_t m = 0; m < task_size + 1; m++) // set boundary condition for t = 0
			{
				double x = (my_rank * common_task_size - 1 + m) * h;
				U_arr_part [0][m] = phi (x);
			}

			for (size_t k = 0; k < K; k++)
			{
				MPI_Recv (&U_arr_part [k + 1][0], 1, MPI_DOUBLE_PRECISION, my_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

				for (size_t m = 1; m < task_size + 1; m++)
				{
					U_arr_part [k + 1][m] = (2 * f ((my_rank * common_task_size + m - 1 - 0.5) * h, (k + 0.5) * tau) + U_arr_part [k][m - 1] * ((1 / tau) + (a / h)) + U_arr_part [k][m] * ((1 / tau) - (a / h)) + U_arr_part [k + 1][m - 1] * ((a / h) - (1 / tau))) / ((1 / tau) + (a / h));
				}

				if ((my_rank * common_task_size + task_size) < M + 1)
				{
					MPI_Send (&U_arr_part [k + 1][task_size], 1, MPI_DOUBLE_PRECISION, my_rank + 1, 0, MPI_COMM_WORLD);
				}
			}

			for (size_t k = 0; k < K + 1; k++) // send parts of solution to the common array
			{
				MPI_Send (&U_arr_part [k][1], task_size, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD);
			}

			//print_array_to_file ("parallel_solution_part.csv", U_arr_part, K + 1, task_size + 1);	
    		free_2D_array (U_arr_part, K + 1);
		}
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

void print_array_to_file (const char * filename, double ** arr, size_t str_num, size_t col_num)
{
	int fd = open (filename, O_WRONLY | O_CREAT | O_TRUNC, 0777); // fd - file descriptor, link to file
														    // O_WRONLY - the file will be opened for writing only
                                                            // O_TRUNC - if the file already exists and is a regular file and the access mode allows writing it will be truncated to length 0
                                                            // O_CREAT - if pathname does not exist, create it as a regular file
    if (fd < 0)
	{
		perror ("Failed to open file for writing");
		exit (EXIT_FAILURE);
	}

	for (size_t i = 0; i < str_num; i++)
	{
		for (size_t j = 0; j < col_num; j++)
		{
			if (dprintf (fd, "%f, ", arr [i][j]) < 0)  // check if something has been written
			{
	        	perror ("Failed to write to file");
	        	close (fd);
	        	exit (EXIT_FAILURE);
    		}
		}
		if (dprintf (fd, "\n") < 0) // and add new line
		{
	       	perror ("Failed to write to file");
	       	close (fd);
	       	exit (EXIT_FAILURE);
    	}
	}

    close (fd);
}