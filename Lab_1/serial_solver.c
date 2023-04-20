#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/unistd.h>
#include <math.h>

// solve differential equation U't + a*U'x = f (x, t)

double psi (double t); // U (0, t) = psi (t)
double phi (double x); // U (x, 0) = phi (t)
double f (double x, double t); // the right part of the differential equation
double U_analytical (double x, double t);
double find_max_error (double ** U_arr, size_t str_num, size_t col_num, double tau, double h);

double ** allocate_2D_array (size_t str_num, size_t col_num);
void free_2D_array (double ** array, size_t str_num);
void print_array_to_file (const char * filename, double ** arr, size_t str_num, size_t col_num);

int main (void)
{
	const double a = 2; // coefficient of the differential equation

	const double T = 1; // maximum value of t axis
	const int K = 20; // maximum number of element of t axis
	const double tau = T / K; // time step

	const double X = 1; // maximum value of x axis
	const int M = 20; // maximum number of element of x axis
	const double h = X / M; // x step

	double ** U_arr = allocate_2D_array (K, M);

	for (int m = 0; m < M + 1; m++) // set boundary condition for t = 0
	{
		double x = m * h;
		U_arr [0][m] = phi (x);
	}

	for (int k = 1; k < K + 1; k++) // set boundary condition for x = 0
	{
		double t = k * tau;
		U_arr [k][0] = psi (t);
	}

	for (int k = 0; k < K; k++) // step-by-step solution
	{
		for (int m = 1; m < M + 1; m++)
		{
			U_arr [k + 1][m] = (2 * f ((m - 0.5) * h, (k + 0.5) * tau) + U_arr [k][m - 1] * ((1 / tau) + (a / h)) + U_arr [k][m] * ((1 / tau) - (a / h)) + U_arr [k + 1][m - 1] * ((a / h) - (1 / tau))) / ((1 / tau) + (a / h));
		}
	}

	printf ("Max error = %E\n", find_max_error (U_arr, K + 1, M + 1, tau, h));
	print_array_to_file ("solution.csv", U_arr, K + 1, M + 1);

    free_2D_array (U_arr, K);
	return 0;
}

double ** allocate_2D_array (size_t str_num, size_t col_num)
{
	double ** arr = (double **) calloc (str_num + 1, sizeof (double *));
    if (arr == NULL)
    {
    	fprintf (stderr, "U_arr pointer is NULL\n");
    	exit (EXIT_FAILURE);
    }
    for (int k = 0; k < str_num + 1; k++)
    {
        arr [k] = calloc (col_num + 1, sizeof (double));
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
	for (int k = 0; k < str_num + 1; k++)
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
		if (dprintf (fd, "\n") < 0) // write the last number in the line and add new line
		{
	       	perror ("Failed to write to file");
	       	close (fd);
	       	exit (EXIT_FAILURE);
    	}
	}

    close (fd);
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