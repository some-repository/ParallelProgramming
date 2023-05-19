#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>

const size_t NUM_THREADS = 10;

typedef struct 
{
	int tid;
} thread_data_t;

void thread_function_helloworld (void * arg)
{
	thread_data_t * thread_data = (thread_data_t *) arg;

	printf ("Hello from thread %d of %ld\n", thread_data -> tid, NUM_THREADS);
	pthread_exit (NULL);
}

int main (void)
{
	pthread_t thread [NUM_THREADS];
	int rc;
	thread_data_t thread_data [NUM_THREADS];

	for (int i = 0; i < NUM_THREADS; i++)
	{
		thread_data [i].tid = i;
		if (rc = pthread_create (&thread [i], NULL, &thread_function_helloworld, &thread_data [i]))
		{
			fprintf (stderr, "Error: pthread_create, rc = %d\n", rc);
			return EXIT_FAILURE;
		}
	}

	for (int i = 0; i < NUM_THREADS; ++i)
	{
		pthread_join (thread [i], NULL);
	}

	return EXIT_SUCCESS;
}