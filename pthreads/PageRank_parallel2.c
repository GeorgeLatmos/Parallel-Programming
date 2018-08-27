#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <pthread.h>

double **alloc_2d_init(int rows, int cols);
double **readDataFromFile();
void PageRank(int lo, int up);
void PageRankParallel(int lo, int up);
void *threads_do(void *thread_data);
double norm(double *x, double *y);

typedef struct 
  {
    int id;
    int lo;
    int up;
   
} THREAD_DATA;

volatile int done;

int N = 4000; 
double **adjacencyMatrix;
double *prev;
double *next;
int maxThreads = 0;
int ActiveThreads = 0;
double tol = 0.001;
double absolute = 0.0;

pthread_barrier_t barr; 

struct timeval startwtime, endwtime;
double seq_time;
double search_time = 0.0;

int main(int argc, char **argv) {

  if (argc < 2) {
    fprintf(stderr, "%s: %s nthread\n", argv[0], argv[0]);
    exit(-1);
  }

  maxThreads = 1<<atoi(argv[1]);

  prev = (double *)malloc(N*sizeof(double));
  next = (double *)malloc(N*sizeof(double));

  adjacencyMatrix = readDataFromFile();

  printf("Starting PageRank Algorithm...\n");

  gettimeofday (&startwtime, NULL);

  PageRankParallel(0,N);

  gettimeofday (&endwtime, NULL);

      seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
        + endwtime.tv_sec - startwtime.tv_sec);

  search_time += seq_time;

  printf("Done in %f seconds!!\n", seq_time);

  return 0;
  
}

/*This is a function that allocates memory in the heap for 2D arrays*/
double **alloc_2d_init(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    int i;
    for (i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

/*This is a function that reads data from a binary file*/
double **readDataFromFile(){

   FILE *fp;
   double *buffer;
   double **ret_buf;
   int size,number_of_elements;
   int i,j,counter;

   fp = fopen("Data.bin","rb");

   fseek(fp, 0, SEEK_END);

   size = ftell(fp); //total size in bytes

   number_of_elements = size/8; //total number of elements

   buffer = (double *)malloc(number_of_elements*sizeof(double));

   ret_buf = alloc_2d_init(N,N);

   fseek(fp, 0, SEEK_SET);

   fread(buffer,sizeof(buffer),number_of_elements,fp);

   fclose(fp);

   counter = 0;
   for(j=0; j<N; j++){
    for(i=0; i<N; i++){
      ret_buf[i][j] = (float)buffer[counter];
      counter++;
    }
   }
   return ret_buf;
}

/*This is the function which is being executed by every thread*/
void *threads_do(void *thread_data){

  THREAD_DATA *tdata = thread_data;

  int tid = (*tdata).id;
  int lo = (*tdata).lo;
  int up = (*tdata).up;

  PageRankParallel(lo,up);
}

/*This is a function that performs PageRank Algorithm using multiple threads*/
void PageRankParallel(int lo, int up){

  int parallel = 0;

  pthread_t threads[2];

  pthread_mutex_t mutexsum;

  pthread_mutex_init(&mutexsum, NULL);

  if(ActiveThreads < maxThreads){

    ActiveThreads += 2;

    pthread_mutex_unlock(&mutexsum);

    THREAD_DATA t1;
    THREAD_DATA t2;
    
    t1.id = 0;
    t1.lo = lo;
    t1.up = lo+up/2;

    t2.id = 1;
    t2.lo = lo+up/2;
    t2.up = up;

    pthread_create(&threads[0], NULL, threads_do, &t1);
    pthread_create(&threads[1], NULL, threads_do, &t2);

    parallel = 1;

  }

  if(parallel){
    
    int i;
        
    for (i=0; i<2; i++){
      pthread_join(threads[i], NULL);
    }

  }else{

    int i;

    PageRank(lo,up);

  }

}

/*This is the PageRank Algorithm*/
void PageRank(int lo, int up){

  int i,j;
  for(i=lo; i<up; i++){
    prev[i] = 0.85+1/N;
  }

  __sync_fetch_and_add(&done, 1);

  for(i=lo; i<up; i++){

      next[i] = 0.0;
      for(j=0; j<N; j++){
         next[i] += adjacencyMatrix[i][j]*prev[j];
      }

    }
}

/*This is a function that calculates the norm of the difference between two vectors*/
double norm(double *x, double *y){
    
  int i;
  double abs = 0;

  for(i=0; i<N; i++){
    abs += pow(x[i]-y[i],2);
  }

  return abs;

}


