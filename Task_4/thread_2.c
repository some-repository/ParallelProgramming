#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>

const size_t NUM_THREADS = 4;

typedef struct 
{
	size_t tid;
	size_t N;
	double * arr_ptr;
} thread_data_t;

void thread_function_sum (void * arg)
{
	thread_data_t * thread_data = (thread_data_t *) arg;
	size_t N = thread_data -> N;
	size_t my_number = (thread_data -> tid);

	size_t epp = 0;
	if (N % NUM_THREADS == 0)
	{
		epp = N / NUM_THREADS; // elements per thread
	}
	else
	{
		epp = (N / NUM_THREADS) + 1; // elements per thread
	}

	double sum = 0;
	for (size_t i = my_number * epp + 1; i < (my_number + 1) * epp + 1; i++)
	{
		if (i > N)
		{
			break;
		}
		sum += 1 / (double) i;
	}
	(thread_data -> arr_ptr) [my_number] = sum;
	pthread_exit (NULL);
}

int main (void)
{
	size_t N = 0;
	scanf ("%ld", &N);
	while (!(N > 0))
	{
		printf ("Wrong N\n");
		scanf ("%ld", &N);
	}

	pthread_t thread [NUM_THREADS];
	int rc;
	thread_data_t thread_data [NUM_THREADS];
	double sum_arr [NUM_THREADS];

	for (size_t i = 0; i < NUM_THREADS; i++)
	{
		thread_data [i].tid = i;
		thread_data [i].arr_ptr = &sum_arr [0];
		thread_data [i].N = N;
		if (rc = pthread_create (&thread [i], NULL, &thread_function_sum, &thread_data [i]))
		{
			fprintf (stderr, "Error: pthread_create, rc = %d\n", rc);
			return EXIT_FAILURE;
		}
	}

	double sum = 0;
	for (size_t i = 0; i < NUM_THREADS; ++i)
	{
		pthread_join (thread [i], NULL);
		sum += sum_arr [i];
		printf ("%f\n", sum_arr [i]);
	}

	printf("N = %ld, sum = %f\n", N, sum);
	return EXIT_SUCCESS;
}