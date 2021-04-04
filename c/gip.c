#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <sys/sysinfo.h>
#include <pthread.h>

#define NSIZE       128
#define VERIFY      0

#define SWAP(a,b)       {double tmp; tmp = a; a = b; b = tmp;}
#define SWAPINT(a,b)       {register int tmp; tmp = a; a = b; b = tmp;}
#define ABS(a)          (((a) > 0) ? (a) : -(a))

void *internal_calc(void *data);

double **matrix,*B,*V,*C;
int *swap;
int num_threads = 1;

/* Allocate the needed arrays */

void allocate_memory(int size)
{
	double *tmp;
	int i;

	matrix = (double**)malloc(size*sizeof(double*));
	assert(matrix != NULL);
	tmp = (double*)malloc(size*size*sizeof(double));
	assert(tmp != NULL);

	for(i = 0; i < size; i++){
		matrix[i] = tmp;
		tmp = tmp + size;
	}

	B = (double*)malloc(size * sizeof(double));
	assert(B != NULL);
	V = (double*)malloc(size * sizeof(double));
	assert(V != NULL);
	C = (double*)malloc(size * sizeof(double));
	assert(C != NULL);
	swap = (int*)malloc(size * sizeof(int));
	assert(swap != NULL);
}

/* Initialize the matirx with some values that we know yield a
 * solution that is easy to verify. A correct solution should yield
 * -0.5, and 0.5 for the first and last C values consecutively, and 0
 * for the rest, though it should work regardless */

void initMatrix(int nsize)
{
	int i,j;
	for(i = 0 ; i < nsize ; i++){
		for(j = 0; j < nsize ; j++) {
			matrix[i][j] = ((j < i )? 2*(j+1) : 2*(i+1));
		}
		B[i] = (double)i;
		swap[i] = i;
	}
}

/* Get the pivot row. If the value in the current pivot position is 0,
 * try to swap with a non-zero row. If that is not possible bail
 * out. Otherwise, make sure the pivot value is 1.0, and return. */

void getPivot(int nsize, int currow)
{
	int i,irow;
	double big;
	double tmp;

	big = matrix[currow][currow];
	irow = currow;

	if (big == 0.0) {
		for(i = currow ; i < nsize; i++){
			tmp = matrix[i][currow];
			if (tmp != 0.0){
				big = tmp;
				irow = i;
				break;
			}
		}
	}

	if (big == 0.0){
		printf("The matrix is singular\n");
		exit(-1);
	}

	if (irow != currow){
		for(i = currow; i < nsize ; i++){
			SWAP(matrix[irow][i],matrix[currow][i]);
		}
		SWAP(B[irow],B[currow]);
		SWAPINT(swap[irow],swap[currow]);
	}


	{
		double pivotVal;
		pivotVal = matrix[currow][currow];

		if (pivotVal != 1.0){
			matrix[currow][currow] = 1.0;
			for(i = currow + 1; i < nsize; i++){
				matrix[currow][i] /= pivotVal;
			}
			B[currow] /= pivotVal;
		}
	}
}


typedef struct thread_data{
	int i;
	int t;
	int pivotVal;
	int nsize;
}thread_t;



/* For all the rows, get the pivot and eliminate all rows and columns
 * for that particular pivot row. */

void computeGauss(int nsize, int num_threads)
{
	int i,j,k;
	double pivotVal;
	int rc;
	void *status;
	pthread_t *threadsInt = (pthread_t*) malloc(num_threads*sizeof(pthread_t));

	for(i = 0; i < nsize; i++){
		getPivot(nsize,i);

		pivotVal = matrix[i][i];
		//spawn threads
		for(int t = 0; t < num_threads; t++){
			thread_t *tmp = malloc(sizeof(thread_t));
			tmp->i = i;
			tmp->t = t;
			tmp->nsize = nsize;
			tmp->pivotVal = pivotVal;
			if(pthread_create(&threadsInt[t],NULL,internal_calc,(void*)tmp)){
				fprintf(stderr,"Failed to create thread\n");
				exit(EXIT_FAILURE);
			}
		}
		//join threads
		for(int t = 0; t<num_threads; t++){
			pthread_join(threadsInt[t], &status);
		}
	}
}
void *internal_calc(void *thread_data){
	thread_t *data = (thread_t *)thread_data;
	int i = data->t;
	int j = data->i;
	int pivotVal = data->pivotVal;
	int nsize = data->nsize;
	int k;
	for (j = i + 1 ; j < nsize; j=j+num_threads){
		pivotVal = matrix[j][i];
		matrix[j][i] = 0.0;
		for (k = i + 1 ; k < nsize; k++){
			matrix[j][k] -= pivotVal * matrix[i][k];
		}
		B[j] -= pivotVal * B[i];
	}
	pthread_exit(NULL);
}


/* Solve the equation. That is for a given A*B = C type of equation,
 * find the values corresponding to the B vector, when B, is all 1's */

void solveGauss(int nsize)
{
	int i,j;

	V[nsize-1] = B[nsize-1];
	for (i = nsize - 2; i >= 0; i --){
		V[i] = B[i];
		for (j = nsize - 1; j > i ; j--){
			V[i] -= matrix[i][j] * V[j];
		}
	}

	for(i = 0; i < nsize; i++){
		C[i] = V[i];//V[swap[i]];
	}
}

extern char * optarg;
int main(int argc,char *argv[])
{
	int i;
	struct timeval start;
	struct timeval finish;
	long compTime;
	double Time;
	int nsize = NSIZE;
	while((i = getopt(argc,argv,"s:n:")) != -1){
		switch(i){
			case 's':
				{
					int s;
					s = atoi(optarg);
					if (s > 0){
						nsize = s;
					} else {
						fprintf(stderr,"Entered size is negative, hence using the default (%d)\n",(int)NSIZE);
					}
				}
				break;
			case 'n':
				{
					int n = atoi(optarg);
					if (n < 1 || n > get_nprocs()){
						fprintf(stderr, "Thats a poor choice of num_threads, try n > 0 or n < %d\n", (int)get_nprocs());
						exit(EXIT_FAILURE);
					}
					else{
						num_threads = n;
					}
				}
				break;
			default:
				assert(0);
				break;
		}
	}

	allocate_memory(nsize);

	gettimeofday(&start, 0);
	initMatrix(nsize);
	computeGauss(nsize,num_threads);
#if VERIFY
	solveGauss(nsize);
#endif
	gettimeofday(&finish, 0);

	compTime = (finish.tv_sec - start.tv_sec) * 1000000;
	compTime = compTime + (finish.tv_usec - start.tv_usec);
	Time = (double)compTime;

	printf("Application time: %f Secs\n",(double)Time/1000000.0);

#if VERIFY
	for(i = 0; i < nsize; i++)
		printf("%6.5f %5.5f\n",B[i],C[i]);
#endif

	return 0;
}
