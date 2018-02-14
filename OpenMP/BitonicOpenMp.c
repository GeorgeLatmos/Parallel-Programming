#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>

struct timeval startwtime, endwtime;
double seq_time;


int N;          // data array size
int *a;         // data array to be sorted


int ActiveThreads = 0;
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
void spawn_threads(int start, int end,int k,int j);
int cmpfuncD(const void * a, const void * b);
int cmpfuncA(const void * a, const void * b);
void recBitonicSortP(int lo, int cnt, int dir);
void bitonicMergeP(int lo, int cnt, int dir);
void threads_merge(int lo,int cnt,int dir);
void spawn_th(int start,int end,int k, int j);
void thread_imp(int start,int end);
void spawn(int start,int end);

/** the main program **/ 
int main(int argc, char **argv) {

  if (argc != 3) {
    printf("Invalid Number of arguments! Aborting...\n");
    exit(1);
  }

  MaxThreads = 1<<atoi(argv[2]);

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

  ActiveThreads = 0;

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

  ActiveThreads = 256;
  
  printf("\nRecursive! Serial...\n");

  init();
  gettimeofday (&startwtime, NULL);
  sort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
          + endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive wall clock time = %f\n\n", seq_time);

  test();

  ActiveThreads = 0;

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

void bitonicMergeParallel(int lo, int cnt, int dir) {

  if (cnt>1) {

    int k = cnt/2;

    omp_lock_t writelock;

    omp_init_lock(&writelock);

    int i;

    for (i=lo; i<lo+k; i++){
        compare(i,i+k,dir);
    }

    threads_merge(lo,cnt,dir);

    omp_set_lock(&writelock);
    ActiveThreads -= 2;
    omp_unset_lock(&writelock);

  }
}


void threads_merge(int lo,int cnt,int dir){

  int k = cnt/2;

  #pragma omp parallel num_threads(2)
  {
        int id = omp_get_thread_num();
        int ilo,icnt;

        if (id == 0){
          ilo = lo;
          icnt = k;
        }else{
          ilo = lo+k;
          icnt = k;
        }

        bitonicMerge(ilo,icnt,dir);

  }
}

/** function recBitonicSort() 
    first produces a bitonic sequence by recursively sorting 
    its two halves in opposite sorting orders, and then
    calls bitonicMerge to make them in the same order 
 **/
void recBitonicSort(int lo, int cnt, int dir) {
  if (cnt>1) {

    if (ActiveThreads<MaxThreads){

      recBitonicSortParallel(lo,cnt,dir);

    }else{
    int k=cnt/2;
    recBitonicSort(lo, k, ASCENDING);
    recBitonicSort(lo+k, k, DESCENDING);
    bitonicMerge(lo, cnt, dir);
    }
  }
}

void recBitonicSortP(int lo, int cnt, int dir) {

  if (cnt>1) {

    if (ActiveThreads<MaxThreads){

      recBitonicSortParallel(lo,cnt,dir);

    }else{

      if (dir == ASCENDING){
        qsort(a+lo,cnt,sizeof(int),cmpfuncA);
    }else{
        qsort(a+lo,cnt,sizeof(int),cmpfuncD);
    }

    }
  }
}

void recBitonicSortParallel(int lo, int cnt, int dir) {
  
  if (cnt>1) {

      int parallel = 0;
      
      omp_lock_t writelock;

      omp_init_lock(&writelock);

      omp_set_lock(&writelock);

      if (ActiveThreads < MaxThreads){

        ActiveThreads += 2;

        omp_unset_lock(&writelock);

        //printf("Active Threads = %d\n",ActiveThreads);

        parallel = 1;
      }
      
      if (parallel){
        
        spawn_RBSP_threads(lo,cnt,dir);

        if (cnt >=8192){

          bitonicMergeParallel(lo,cnt,dir);

        }else{

          bitonicMerge(lo,cnt,dir);

        }
      

      }else{

        recBitonicSortP(lo,cnt,dir);
        
      }

  }
}

void spawn_RBSP_threads(int lo, int cnt, int dir){

  #pragma omp parallel num_threads(2)
  {

    int k = cnt/2;

    int id = omp_get_thread_num();
    
    if (id == 0){
      recBitonicSortParallel(lo,k,ASCENDING);
      
    }else{
      recBitonicSortParallel(lo+k,k,DESCENDING);
      
    }
  }

}

/** function sort() 
   Caller of recBitonicSort for sorting the entire array of length N 
   in ASCENDING order
 **/
void sort() {
  recBitonicSort(0, N, ASCENDING);
}

void sortParallel() {
  omp_set_nested(1);
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

  omp_set_nested(1);

  spawn(0,N);
}

void spawn(int start,int end){

  int i,j,k;

  int parallel = 0;

  omp_lock_t writelock;

  omp_init_lock(&writelock);

  omp_set_lock(&writelock);

  if (ActiveThreads<MaxThreads){

    ActiveThreads += 2;

    omp_unset_lock(&writelock);

    parallel = 1;
  }

  if (parallel){

    thread_imp(start,end);

    k = (end-start);

    for (j=k>>1; j>0; j=j>>1) {
      spawn_th(start,end-start,k,j);
      
    }

    omp_set_lock(&writelock);
    ActiveThreads -= 2;
    omp_unset_lock(&writelock);

  }else{

  for (k=2; k<=(end-start); k=2*k) {
    for (j=k>>1; j>0; j=j>>1) {
      for (i=start; i<(end-start); i++) {
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
 
}

void thread_imp(int start,int end){

  #pragma omp parallel num_threads(2)
  {
    int id = omp_get_thread_num();
    int istart,iend;

    if (id == 0){
      istart = start;
      iend = start+(end-start)/2;
      spawn(istart,iend);
    }else{
      istart = start+(end-start)/2;
      iend = end;
      qsort(a+istart,iend-istart,sizeof(int),cmpfuncD);
    }

  }

}

void spawn_th(int start,int end,int k, int j){

  #pragma omp parallel num_threads(2)
  {
    int id = omp_get_thread_num();
    int istart,iend;
    int i;

    if (id == 0){
      istart = start;
      iend = start+(end-start)/2;
    }else{
      istart = start+(end-start)/2;
      iend = end;
    }

    for (i=istart; i<iend; i++) {
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







