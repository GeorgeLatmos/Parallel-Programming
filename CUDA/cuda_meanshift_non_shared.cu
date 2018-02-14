#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

const int N = 600;
const int D = 2;
const int blocksize = 4;
const int gridsize = 4;

struct timeval startwtime, endwtime;
double seq_time;

double **alloc_2d_init(int rows, int cols);
double **init();

__device__
void findNext(double *a,double *b, double *c,int index){
	
	double numinator;
    double denuminator;
    double distance;
    double sigma = 1.0;
    double c1;
    int j,k,d;
	
    for(j=0; j<D; j++){
       	numinator = 0.0;
       	denuminator = 0.0;
       	for(k=0; k<N; k++){
           distance = 0.0;
    	   for(d=0; d<D; d++){
    		  distance += (b[index*D+d]-a[k*D+d])*(b[index*D+d]-a[k*D+d]);
    	   }
         if(distance<sigma*sigma){
           c1 = exp(-distance/(2*sigma*sigma));
           denuminator += c1;
           numinator += c1*a[k*D+j];
         }
    	   
    	}
    	c[index*D+j] = numinator/denuminator;
    }
}

__global__
void meanshiftNonShared(double *a, double *b, double *c){

	int index = blockIdx.x*blockDim.x + threadIdx.x;
	int i,j;

  double abs_m;

  double epsilon = 0.0001;

  for(i=index*N/(blockDim.x*blockDim.x); i<(index+1)*N/(blockDim.x*blockDim.x); i++){
        
    abs_m = 1.0;

    while(abs_m > epsilon){

    	findNext(a,b,c,i);

    	abs_m = 0.0;
    	for(j=0; j<D; j++){
    	   abs_m += (c[i*D+j]-b[i*D+j])*(c[i*D+j]-b[i*D+j]);
      }
      for(j=0; j<D; j++){
          b[i*D+j] = c[i*D+j];
      }

    }
  }
}

int main(){

	double **points;
	double **new_points;

	double *lin_points, *lin_prev, *lin_next;

	lin_points = (double *)malloc(N*D*sizeof(double));
	lin_prev = (double *)malloc(N*D*sizeof(double));
	lin_next = (double *)malloc(N*D*sizeof(double));

	new_points = alloc_2d_init(N,D);

	points = init();

	int i,j,index;

	for(i=0; i<N; i++){
		for(j=0; j<D; j++){
			index = j+i*D;
			lin_points[index] = points[i][j];
			lin_prev[index] = points[i][j];
		}
	}

	double *ad, *bd, *cd;
	const int size = N*D*sizeof(double);

	cudaMalloc((void **)&ad,size);
	cudaMalloc((void **)&bd,size);
	cudaMalloc((void **)&cd,size);

	cudaMemcpy(ad,lin_points,size,cudaMemcpyHostToDevice);
	cudaMemcpy(bd,lin_prev,size,cudaMemcpyHostToDevice);

	dim3 dimBlock(blocksize);
	dim3 dimGrid(gridsize);
    
    gettimeofday (&startwtime, NULL);

	meanshiftNonShared<<<dimGrid,dimBlock>>>(ad,bd,cd);

	cudaDeviceSynchronize();

	gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
        + endwtime.tv_sec - startwtime.tv_sec);

  printf("Meanshift Non Shared wall clock time = %f\n\n", seq_time);

  cudaMemcpy(lin_next,cd,size,cudaMemcpyDeviceToHost);
  
	for(i=0; i<N; i++){
		for(j=0; j<D; j++){
			new_points[i][j] = lin_next[j+i*D];
		}
	}

  cudaFree(ad); 
  cudaFree(bd);
  cudaFree(cd);

	return EXIT_SUCCESS;
}

double **init(){

   FILE *fp;
   double *buffer;
   double **ret_buf;
   int size,number_of_elements;
   int i,j,counter;

   fp = fopen("data.bin","rb");

   fseek(fp, 0, SEEK_END);

   size = ftell(fp); //total size in bytes

   number_of_elements = size/8; //total number of elements

   buffer = (double *)malloc(number_of_elements*sizeof(double));

   ret_buf = alloc_2d_init(N,D);

   fseek(fp, 0, SEEK_SET);

   fread(buffer,sizeof(buffer),number_of_elements,fp);

   fclose(fp);

   counter = 0;
   for(j=0; j<D; j++){
    for(i=0; i<N; i++){
      ret_buf[i][j] = buffer[counter];
      counter++;
    }
   }
   return ret_buf;
}

double **alloc_2d_init(int rows, int cols){
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    int i;
    for (i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}