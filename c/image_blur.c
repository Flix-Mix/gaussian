// The program takes a color image as input and computes the greyscale intensity for each pixel.
// Each pixel is classified as dark/bright/none based on its intensity. 
// Each pixel is then blurred with its surroundings.  
// Created by: Sharanyan Srikanthan, 1/19/2016
// Copyright: University of Rochester
#include <pthread.h>
#include <stdio.h>
#include <sched.h>
#include <stdlib.h>
#include <time.h>
#include "hrtimer_x86.h"

//Defining constants
#define MAX_VAL 254
#define FILTER_LEN 2
#define INTENSITY_THRESHOLD 10
#define ZERO_PROB 25

const float r_w = 0.2989;
const float g_w = 0.5870;
const float b_w = 0.1140;
int num_threads = 32;

extern char *optarg;
// Arrays for storing matrices
int m = 1; 
int **r, **g, **b; 
float **intensity, **res;
// Per-thread counters
int *zero_count;
int *intensity_count;
pthread_t *threadsInt,*threadsRes;

// Deallocates matrices
void dealloc()
{
	int len;
	for (len = 0; len < m; len++)
	{
		free(r[len]);
		free(g[len]);
		free(b[len]);
		free(res[len]);
		free(intensity[len]);
	}
	free(r);free(g);free(b);free(res);free(intensity);
	free(threadsInt);free(threadsRes);
	free(zero_count);free(intensity_count);
	return;
}

// Calculate intensity of a given element using R,G,B components
// Counts dark (r = g = b = 0) and bright (intensity > threshold defined) elements
void *IntensityCalc(void *threadid)
{
	long tid;
	tid = (long)threadid;
	int i, j;

	// Each thread is responsible for one row at a time. 
	for (i = tid; i < m; i=i+num_threads)
		for (j = 0; j < m; j++)
		{
			intensity[j][i] = (r_w * r[j][i]) + (g_w * g[j][i]) + (b_w * b[j][i]);
			if (r[j][i] == 0)
				zero_count[tid]++;
			if (intensity[j][i] > INTENSITY_THRESHOLD)
				intensity_count[tid]++;
		}
	pthread_exit(NULL);
}
void VerifyIntensityCalc()
{
	int i, j;
	float temp;

	for (i = 0; i < m; i++)
		for (j = 0; j < m; j++)
		{
			temp = (r_w * r[j][i]) + (g_w * g[j][i]) + (b_w * b[j][i]);
			if (intensity[j][i] != temp)
			{
				printf("Error in intensity calculation at %d , %d \n",j,i);
				return;
			}
		}
	return;
}
void *FilterCalc(void *threadid)
{
	long tid;
	tid = (long)threadid;
	int i, j, len, bth; // Loop counter variables
	int left, right, up, down;

	// Each thread is responsible for one column at a time.
	for (i = tid; i < m; i=i+num_threads)
	{
		for (j = 0; j < m; j++)
		{
			res[j][i] = 0;
			// Determine boundaries to apply filter to avoid segmentation fault
			left = ((i - FILTER_LEN) > 0)?(i - FILTER_LEN):0;
			right = ((i + FILTER_LEN) < (m-1))?(i + FILTER_LEN):(m-1);
			up = ((j - FILTER_LEN) > 0)?(j - FILTER_LEN):0;
			down = ((j + FILTER_LEN) < (m-1))?(j + FILTER_LEN):(m-1);
			
			// Perform accumulation and averaging
			for (len = up; len <= down;  len++)
				for (bth = left; bth <= right; bth++)
					res[j][i] += intensity[bth][len];
			res[j][i] = res[j][i]/((right-left+1)*(down-up+1));
		}
	}
	pthread_exit(NULL);
}
void VerifyFilterCalc()
{
	int i, j, len, bth; // Loop counter variables
	int left, right, up, down;
	float temp;

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			temp = 0;
			// Determine boundaries to apply filter to avoid segmentation fault
			left = ((i - FILTER_LEN) > 0)?(i - FILTER_LEN):0;
			right = ((i + FILTER_LEN) < (m-1))?(i + FILTER_LEN):(m-1);
			up = ((j - FILTER_LEN) > 0)?(j - FILTER_LEN):0;
			down = ((j + FILTER_LEN) < (m-1))?(j + FILTER_LEN):(m-1);
			
			// Perform accumulation and averaging
			for (len = up; len <= down;  len++)
				for (bth = left; bth <= right; bth++)
					temp += intensity[bth][len];
			temp = temp/((right-left+1)*(down-up+1));
			if (temp != res[j][i])
			{
				printf("Error in filtering at %d, %d \n",j,i);
				return;
			}
		}
	}
	return;
}
int main (int argc, char *argv[])
{
	int rc, verify;
	long t;
	void *status;
	double start_time, mid_time, end_time;
	int len, bth;
	int zero_c = 0, intensity_c = 0;

	while ((rc = getopt (argc, argv, "vhm:n:")) != -1)
		switch (rc) 
		{
			case 'm':
				m = atoi (optarg);
				break;
			case 'n':
				num_threads = atoi (optarg);
				break;
			case 'v':
				verify = 1;
				break;
			case 'h':
				printf("Usage: ./program -n <num threads> -m <matrix size>\n");
				return 0;
				break;
		}

	printf("\n Matrix Size %d Threads %d Verification %d\n", m,num_threads,verify);

	// Matrix allocation and initialization
	srand(time(NULL));
	r = (int**) malloc(m*sizeof(int*));
	g = (int**) malloc(m*sizeof(int*));
	b = (int**) malloc(m*sizeof(int*));
	intensity = (float**) malloc(m*sizeof(float*));
	res = (float**) malloc(m*sizeof(float*));
	intensity_count = (int*) malloc(num_threads * sizeof(int));
	zero_count = (int*) malloc(num_threads * sizeof(int));
	threadsInt = (pthread_t*) malloc(num_threads * sizeof(pthread_t));
	threadsRes = (pthread_t*) malloc(num_threads * sizeof(pthread_t));

	for (len = 0; len < m; len++)
	{
		r[len] = (int*) malloc(m*sizeof(int));
		g[len] = (int*) malloc(m*sizeof(int));
		b[len] = (int*) malloc(m*sizeof(int));
		intensity[len] = (float*) malloc(m*sizeof(float));
		res[len] = (float*) malloc(m*sizeof(float));
		for (bth = 0; bth < m; bth++)
		{
			if (rand() % ZERO_PROB == 0)
			{
				r[len][bth] = 0;
				g[len][bth] = 0;
				b[len][bth] = 0;
			}
			else
			{
				r[len][bth] = 1 + (rand() % MAX_VAL);
				g[len][bth] = 1 + (rand() % MAX_VAL);
				b[len][bth] = 1 + (rand() % MAX_VAL);
			}
			res[len][bth] = 0;
			intensity[len][bth] = 0;
		}
	}
	printf("Matrix has been initialized \n");
	start_time = gethrtime_x86();
	// Spawn threads for intensity calculation
	for(t=0; t<num_threads; t++)
	{
		rc = pthread_create(&threadsInt[t], NULL, IntensityCalc, (void *)(t));
		if (rc)
		{
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			dealloc();
			return -1;
		}
	}	 
	// Wait for thread termination
	for(t=0; t<num_threads; t++)
		pthread_join(threadsInt[t], &status);
	for (t = 0; t < num_threads; t++)
	{
		zero_c += zero_count[t];
		intensity_c += intensity_count[t];
	}
	
	mid_time = gethrtime_x86();
	printf("Intensity done in %f secs \nIntensity Count: %d Zero Count: %d \n",mid_time - start_time, intensity_c,zero_c);
	// Spawn threads for filtering
	for(t=0; t<num_threads; t++)
	{
		rc = pthread_create(&threadsRes[t], NULL, FilterCalc, (void *)(t));
		if (rc)
		{
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			dealloc();
			return -1;
		}
	}
	   
	// Wait for thread termination
	for(t=0; t<num_threads; t++)
		pthread_join(threadsRes[t], &status);
	end_time = gethrtime_x86();
	printf("Filtering done in %f \n",end_time - mid_time);
	if (verify == 1)
	{
		printf("Starting intensity verification \n");
		VerifyIntensityCalc();
		printf("Intensity verification is complete\n");
		printf("Starting filter verification \n");
		VerifyFilterCalc();
		printf("Filter verification is complete \n");
	}
		

	dealloc();
	return 0;
}
