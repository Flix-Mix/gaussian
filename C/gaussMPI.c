#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <mpi.h>
#include <sys/sysinfo.h>
#include <math.h>

#define NSIZE       128
#define VERIFY      0

#define SWAP(a,b)       {double tmp; tmp = a; a = b; b = tmp;}
#define SWAPINT(a,b)       {register int tmp; tmp = a; a = b; b = tmp;}
#define ABS(a)          (((a) > 0) ? (a) : -(a))

double **matrix,*B,*V,*C;
double **matbuf; 
double *Bbuf;//,*bp,*bpp;
int *swap;
int num_procs=1;
int id,num_processors;
int *sendcounts,*bsendcounts;
int *displ,*displ2;

//void fill(int v, int i, int nsize, int fg);
void fill_sendcounts(int nsize,int i);
void fill_displacement(int nsize,int i);


/* Allocate the needed arrays */

void allocate_memory(int size)
{
	double *tmp,*tmp2;
	int i;

	matbuf = (double**)malloc(size*sizeof(double*));

	matrix = (double**)malloc(size*sizeof(double*));
	assert(matrix != NULL);
	tmp = (double*)malloc(size*size*sizeof(double));
	tmp2 = (double*)malloc(size*size*sizeof(double));
	assert(tmp != NULL);

	for(i = 0; i < size; i++){
		matrix[i] = tmp;
		tmp = tmp + size;

		matbuf[i] = tmp2;
		tmp2 = tmp2 + size;
	}

	B =    (double*)malloc(size*sizeof(double));
	Bbuf = (double*)malloc(size*sizeof(double));
//	bp = (double*)malloc(size*sizeof(double));
//	bpp = (double*)malloc(size*sizeof(double));

	displ = malloc(size*sizeof(int));
	displ2 = malloc(size*sizeof(int));
	sendcounts = malloc(size*sizeof(int));
	bsendcounts = malloc(size*sizeof(int));


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


void fill_displacement(int nsize,int i){
	for(int l = 0; l<num_processors; l++){
		if(l != 0) { displ[l] = displ[l-1] + sendcounts[l-1]; }
		else{displ[l] = 0;}

		if(l != 0) { displ2[l] = displ2[l-1] + bsendcounts[l-1]; }
		else{displ2[l] = 0;}

		//printf("displ: %d\n", displ[l]);
		//printf("displ2: %d\n", displ2[l]);
	}
//	printf("matrix displ: [");
//	for(int j = 0; j<num_processors; j++){
//		printf(" %d,", displ[j]);
//	}
//	printf("]\nB displ: [");
//	for(int j = 0; j<num_processors; j++){
//		printf(" %d,", displ2[j]);
//	}
//	printf("]\n");

}

/* For all the rows, get the pivot and eliminate all rows and columns
 * for that particular pivot row. */

void computeGauss(int nsize)
{
	int i,j,k;
	double pivotVal;
	for(i = 0; i < nsize; i++){
		getPivot(nsize,i);

		pivotVal = matrix[i][i];
		/* MPI */
		MPI_Bcast(matrix[i],nsize,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&B[i],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		fill_sendcounts(nsize,i);
		fill_displacement(nsize,i);
		if(MPI_Scatterv((const void *)matrix[i+1],sendcounts,displ,MPI_DOUBLE,matbuf[i+1],sendcounts[id],MPI_DOUBLE,0,MPI_COMM_WORLD)  != MPI_SUCCESS){
			fprintf(stderr,"%s\n","Failed to scatter");
			exit(EXIT_FAILURE);
		}
		if(MPI_Scatterv(&B[i+1],bsendcounts,displ2,MPI_DOUBLE,&Bbuf[i+1],bsendcounts[id],MPI_DOUBLE,0,MPI_COMM_WORLD)  != MPI_SUCCESS){
			fprintf(stderr,"%s\n","Failed to scatter");
			exit(EXIT_FAILURE);
		}
		/* MPI */
		for (j = i + 1 ; j <= bsendcounts[id]; j++){
			pivotVal = matbuf[j][i];
			matbuf[j][i] = 0.0;
			for (k = i + 1 ; k < nsize; k++){
				matbuf[j][k] -= pivotVal * matrix[i][k];
			}
			Bbuf[j] -= pivotVal * B[i];
		}

		/* MPI */
		if(MPI_Gatherv(matbuf[i+1],sendcounts[id],MPI_DOUBLE,matrix[i+1],sendcounts,displ, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS){
			fprintf(stderr,"%s\n","Failed to gather");
			exit(EXIT_FAILURE);
		}
		if(MPI_Gatherv(&Bbuf[i+1],bsendcounts[id],MPI_DOUBLE,&B[i+1],bsendcounts,displ2, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS){
			fprintf(stderr,"%s\n","Failed to gather");
			exit(EXIT_FAILURE);
		}
		/* MPI */
	}
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

void fill_sendcounts(int nsize,int i){
	int rtd = nsize-i-1;
	int rltd = rtd;
	int nrpp = ceil((double)rtd/num_processors);
	for(int n = 0; n<num_processors; n++){
		if(rltd >= nrpp){
			sendcounts[n] = nrpp*nsize;
			bsendcounts[n] = nrpp;
			rltd -= nrpp;
		}
		else if(rltd == 0){
			sendcounts[n] = 0;
			bsendcounts[n] = 0;
		}
		else{
			sendcounts[n] = rltd*nsize;
			bsendcounts[n] = rltd;
			rltd -= rltd;
		}
	}

//	for(int i = 0; i<num_processors; i++){
//		printf("sc: %d\n", sendcounts[i]);
//		printf("bc: %d\n", bsendcounts[i]);
//	}
}

extern char * optarg;
int main(int argc,char *argv[])
{
	MPI_Init(&argc,&argv);
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
					int n;
					n = atoi(optarg);
					if (n<1 || n> get_nprocs()){
						fprintf(stderr, "Poor choice of num threads, try n > 0 or n < %d\n", (int)get_nprocs());
						exit(EXIT_FAILURE);
					}
					else{
						num_procs = n;
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






	/*
	MPI
	*/
//	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
	//printf("BEFORE FILL\n");
	//fill_sendcounts(nsize);
	//printf("AFTER FILL\n");
//	for(int i = 0; i<num_processors; i++){
//		printf("%d\n",sendcounts[i]);
//	}
//	printf("DONE\n");
	computeGauss(nsize);

	MPI_Finalize();
	/*
	MPI
	*/






#if VERIFY
	solveGauss(nsize);
#endif
	gettimeofday(&finish, 0);

	compTime = (finish.tv_sec - start.tv_sec) * 1000000;
	compTime = compTime + (finish.tv_usec - start.tv_usec);
	Time = (double)compTime;

	//printf("Application time: %f Secs\n",(double)Time/1000000.0);
	printf("%f\n",(double)Time/1000000.0);

#if VERIFY
	for(i = 0; i < nsize; i++)
		printf("%6.5f %5.5f\n",B[i],C[i]);
#endif

	return 0;
}
