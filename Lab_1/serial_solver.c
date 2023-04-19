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

void print_array_to_file (const char * filename, double ** arr, size_t str_num, size_t col_num);

int main (void)
{
	const double a = 2; // coefficient of the differential equation

	const double T = 1; // maximum value of t axis
	const int K = 10 * 25; // maximum number of element of t axis
	const double tau = T / K; // time step

	const double X = 1; // maximum value of x axis
	const int M = 10 * 10; // maximum number of element of x axis
	const double h = X / M; // x step

	double ** U_arr = (double **) calloc (K + 1, sizeof (double *));
    if (U_arr == NULL)
    {
    	fprintf (stderr, "U_arr pointer is NULL\n");
    	return 1;
    }
    for (int k = 0; k < K + 1; k++)
    {
        U_arr [k] = calloc (M + 1, sizeof (double));
        if (U_arr [k] == NULL)
    	{
    		fprintf (stderr, "U_arr pointer is NULL\n");
    		return 1;
    	}
    }

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
			U_arr [k + 1][m] = 2 * f ((m - (1 / 2)) * h, (k + (1 / 2)) * tau) + U_arr [k][m - 1] * ((1 / tau) + (a / h)) + U_arr [k][m] * ((1 / tau) - (a / h)) + U_arr [k + 1][m - 1] * ((a / h) - (1 / tau));
		}
	}

	print_array_to_file ("solution.txt", U_arr, K + 1, M + 1);

	for (int k = 0; k < K + 1; k++)
    {
        free (U_arr [k]);
    }
    free (U_arr);

	return 0;
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
			if (dprintf (fd, "%E ", arr [i][j]) < 0)  // check if something has been written
			{
	        	perror ("Failed to write to file");
	        	close (fd);
	        	exit (EXIT_FAILURE);
    		}
		}
		if (dprintf (fd, "\n") < 0) // add new line
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
	return 0;
}