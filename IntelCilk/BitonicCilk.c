#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>


struct timeval startwtime, endwtime;
double seq_time;


int N;          // data array size
int *a;         // data array to be sorted

//int ActiveThreads = 0;

int MaxThreads = 0;

const int ASCENDING  = 1;
const int DESCENDING = 0; 


void init(void);
void print(void);
void sort(void);
void sortParallel(void);
void test(void);
void exchange1(int i, int j);
void compare(int i, int j, int dir);
void bitonicMerge(int lo, int cnt, int dir);
void recBitonicSort(int lo, int cnt, int dir);
void recBitonicSortParallel(int lo, int cnt, int dir);
void impBitonicSort(void);
void impBitonicSortParallel();
void spawn_RBSP_threads();
int cmpfuncA(const void * a, const void * b);
int cmpfuncD(const void * a, const void * b);
void spawn_threads(int start, int end,int k,int j);
void bitonicMergeParallel(int lo, int cnt, int dir);
void thread_imp_for(int start,int end,int k,int j);
void spawn(int start,int end);
void sort_imp(int start,int end);

/** the main program **/ 
int main(int argc, char **argv) {

  if (argc != 2) {
    printf("Invalid Number of arguments! Aborting...\n");
    exit(1);
  }

  //MaxThreads = 1<<atoi(argv[2]);

  N = 1<<atoi(argv[1]);
  a = (int *) malloc(N * sizeof(int));

  printf("Imperative! Serial...\n");
  
  init();

  gettimeofday (&startwtime, NULL);
  impBitonicSort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Imperative wall clock time = %f\n\n", seq_time);

  test();
  //print();

  printf("\nImperative! Parallel...\n");

  init();

  gettimeofday (&startwtime, NULL);
  impBitonicSortParallel();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Imperative (Parallel) wall clock time = %f\n\n", seq_time);

  //print();
  test();
  
  printf("\nRecursive! Serial...\n");

  init();
  gettimeofday (&startwtime, NULL);
  sort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
          + endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive wall clock time = %f\n\n", seq_time);

  test();

  printf("\nRecursive! Parallel...\n");

  init();
  gettimeofday (&startwtime, NULL);
  sortParallel();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
          + endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive (Parallel) wall clock time = %f\n\n", seq_time);

  test();

  printf("\nqsort!...\n");

  init();
  gettimeofday (&startwtime, NULL);
  qsort(a,N,sizeof(int),cmpfuncA);
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
          + endwtime.tv_sec - startwtime.tv_sec);

  printf("Qsort wall clock time = %f\n\n", seq_time);

  test();
  
}

/** -------------- SUB-PROCEDURES  ----------------- **/ 

int cmpfuncD(const void * a, const void * b) {
   return ( -*(int*)a + *(int*)b );
}

int cmpfuncA(const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

/** procedure test() : verify sort results **/
void test() {
  int pass = 1;
  int i;
  for (i = 1; i < N; i++) {
    pass &= (a[i-1] <= a[i]);
  }

  printf(" TEST %s\n",(pass) ? "PASSed" : "FAILed");
}


/** procedure init() : initialize array "a" with data **/
void init() {
  int i;
  for (i = 0; i < N; i++) {
    a[i] = rand() % N; // (N - i);
  }
}

/** procedure  print() : print array elements **/
void print() {
  int i;
  for (i = 0; i < N; i++) {
    printf("%d\n", a[i]);
  }
  printf("\n");
}


/** INLINE procedure exchange() : pair swap **/
void exchange1(int i, int j) {
  int t;
  t = a[i];
  a[i] = a[j];
  a[j] = t;
}

/** procedure compare() 
   The parameter dir indicates the sorting direction, ASCENDING 
   or DESCENDING; if (a[i] > a[j]) agrees with the direction, 
   then a[i] and a[j] are interchanged.
**/
void compare(int i, int j, int dir) {
  if (dir==(a[i]>a[j])) 
    exchange1(i,j);
}

/** Procedure bitonicMerge() 
   It recursively sorts a bitonic sequence in ascending order, 
   if dir = ASCENDING, and in descending order otherwise. 
   The sequence to be sorted starts at index position lo,
   the parameter cbt is the number of elements to be sorted. 
 **/

void bitonicMerge(int lo, int cnt, int dir) {
  if (cnt>1) {
    int k=cnt/2;
    int i;

    for (i=lo; i<lo+k; i++) 
      compare(i, i+k, dir);
    bitonicMerge(lo, k, dir);
    bitonicMerge(lo+k, k, dir);
  }
}

void bitonicMergeParallel(int lo, int cnt, int dir){
  
  int k=cnt/2;

  int i;

  for (i=lo; i<lo+k; i++){
     compare(i,i+k,dir);
  }

  cilk_spawn bitonicMerge(lo,k,dir);
  bitonicMerge(lo+k,k,dir);
  cilk_sync;



} 


/** function recBitonicSort() 
    first produces a bitonic sequence by recursively sorting 
    its two halves in opposite sorting orders, and then
    calls bitonicMerge to make them in the same order 
 **/
void recBitonicSort(int lo, int cnt, int dir) {
  if (cnt>1) {
    int k=cnt/2;
    recBitonicSort(lo, k, ASCENDING);
    recBitonicSort(lo+k, k, DESCENDING);
    bitonicMerge(lo, cnt, dir);
  }
}

void recBitonicSortP(int lo, int cnt, int dir) {
  if (cnt>1) {

    if (dir == ASCENDING){
        qsort(a+lo,cnt,sizeof(int),cmpfuncA);
    }else{
        qsort(a+lo,cnt,sizeof(int),cmpfuncD);
    }
    

  }
}
void recBitonicSortParallel(int lo, int cnt, int dir) {
  
  if (cnt>1) {

      int k = cnt/2;

      cilk_spawn recBitonicSortP(lo,k,ASCENDING);
      recBitonicSortP(lo+k,k,DESCENDING);
      cilk_sync;

      bitonicMergeParallel(lo,cnt,dir);

  }
}

void spawn_RBSP_threads(int lo, int cnt, int dir){
   
  int k = cnt/2;

  cilk_spawn recBitonicSortParallel(lo,k,ASCENDING);
  recBitonicSortParallel(lo+k,k,DESCENDING);
  cilk_sync;

  bitonicMerge(lo,cnt,dir); 

}

/** function sort() 
   Caller of recBitonicSort for sorting the entire array of length N 
   in ASCENDING order
 **/
void sort() {
  recBitonicSort(0, N, ASCENDING);
}

void sortParallel() {
  
  recBitonicSortParallel(0, N, ASCENDING);
}

/*
  imperative version of bitonic sort
*/
void impBitonicSort() {

  int i,j,k;
  
  for (k=2; k<=N; k=2*k) {
    for (j=k>>1; j>0; j=j>>1) {
      for (i=0; i<N; i++) {
  int ij=i^j;
  if ((ij)>i) {
    if ((i&k)==0 && a[i] > a[ij]) 
        exchange1(i,ij);
    if ((i&k)!=0 && a[i] < a[ij])
        exchange1(i,ij);
  }
      }
    }
  }
}

void impBitonicSortParallel() {

  spawn(0,N);

}

void spawn(int start,int end){

  int i,j,k;

  int istart = start+(end-start)/2;
  int iend = start+(end-start)/2;

  cilk_spawn sort_imp(start,iend);
  qsort(a+istart,end-istart,sizeof(int),cmpfuncD);
  cilk_sync;

  k = end-start;

  for (j=k>>1; j>0; j=j>>1) {
    spawn_threads(start,end-start,k,j);  
  }
 
}

void sort_imp(int start,int end){

  int i,j,k;

  for (k=2; k<=(end-start); k=2*k) {
    for (j=k>>1; j>0; j=j>>1) {

      cilk_for (i=start; i<(end-start); i++) {
  int ij=i^j;
  if ((ij)>i) {
    if ((i&k)==0 && a[i] > a[ij]) 
        exchange1(i,ij);
    if ((i&k)!=0 && a[i] < a[ij])
        exchange1(i,ij);
  }
      }
    
    }
  }

}

void spawn_threads(int start, int end,int k,int j){

  int istart = start+(end-start)/2;
  int iend = start+(end-start)/2;

  cilk_spawn thread_imp_for(start,iend,k,j);
  thread_imp_for(istart,end,k,j);
  cilk_sync;
}

void thread_imp_for(int start,int end,int k,int j){

  int i;

  cilk_for (i=start; i<end; i++) {
    int ij=i^j;
    if ((ij)>i) {
      if ((i&k)==0 && a[i] > a[ij]) 
        exchange1(i,ij);
      if ((i&k)!=0 && a[i] < a[ij])
        exchange1(i,ij);
    }
  }

}