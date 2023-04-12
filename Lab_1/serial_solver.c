#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// solve differential equation U't + a*U'x = f (x, t)

double psi (double t); // U (0, t) = psi (t)
double phi (double x); // U (x, 0) = phi (t)
double f (double x, double t); // the right part of differential equation

double U_analytical (double x, double t);

int main (void)
{
	const double T = 10; // maximum value of t axis
	const int K = 10 * 1000; // maximum number of element of t axis
	const double tau = T / K; // time step

	const double X = 100; // maximum value of x axis
	const int M = 10 * 1000; // maximum number of element of x axis
	const double h = X / M; // x step

	const double a = 0.1; // coefficient of differential equation

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

	for (int k = 0; k < K + 1; k++)
    {
        free (U_arr [k]);
    }
    free (U_arr);

	return 0;
}

double psi (double t)
{
	return 2 * t;
}

double phi (double x)
{
	return 3 * x;
}

double f (double x, double t)
{
	return 0;
}

double U_analytical (double x, double t)
{

}