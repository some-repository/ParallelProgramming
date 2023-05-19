#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>

const size_t NUM_THREADS = 4;

pthread_mutex_t lock_x;
int shared_var = 0;

typedef struct 
{
	size_t tid;
} thread_data_t;

void thread_function (void * arg)
{
	thread_data_t * thread_data = (thread_data_t *) arg;

	printf ("My thread number = %ld of %ld\n", thread_data -> tid, NUM_THREADS);

	pthread_mutex_lock (&lock_x);
	shared_var ++;
	printf ("shared_var = %d\n", shared_var);
	pthread_mutex_unlock (&lock_x);

	pthread_exit (NULL);
}

int main (void)
{
	pthread_t thread [NUM_THREADS];
	int rc;
	thread_data_t thread_data [NUM_THREADS];
	shared_var = 1;

	pthread_mutex_init (&lock_x, NULL);

	for (int i = 0; i < NUM_THREADS; i++)
	{
		thread_data [i].tid = i;
		if (rc = pthread_create (&thread [i], NULL, &thread_function, &thread_data [i]))
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