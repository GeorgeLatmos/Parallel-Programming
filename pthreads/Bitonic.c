#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

typedef struct 
  {
    int id;
    int lo;
    int cnt;
    int dir;  
  } THREAD_DATA;

typedef struct 
  {
    int id;
    int start;
    int end;
    int k;
    int j;  
  } THREAD_DATA_IMP;

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
void recBitonicSortP(int lo, int cnt, int dir);
void impBitonicSort(void);
void impBitonicSortParallel();
void *threads_do(void *threadid);
int cmpfuncD(const void * a, const void * b);
int cmpfuncA(const void * a, const void * b);
void *threads_merge(void *thread_data);
void bitonicMergeParallel(int lo, int cnt, int dir);
void spawn(int start, int end);
void *thread_imp(void *t_data);
void spawn_th(int start,int end,int k, int j);
void *thread_imp_for(void *t_data);


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
  //print();

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

  test();
  //print();
  
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

    pthread_t threads[2];

    pthread_mutex_t mutexsum;

    pthread_mutex_init(&mutexsum, NULL);

    int i;

    for (i=lo; i<lo+k; i++){
        compare(i,i+k,dir);
    }

    THREAD_DATA t1;
    THREAD_DATA t2;
        
    t1.id = 0;
    t1.lo = lo;
    t1.cnt = k;
    t1.dir = dir;
  
    t2.id = 1;
    t2.lo = lo+k;
    t2.cnt = k;
    t2.dir = dir;

    pthread_create(&threads[0], NULL, threads_merge, &t1);
    pthread_create(&threads[1], NULL, threads_merge, &t2);
        
    for (i=0; i<2; i++){
      pthread_join(threads[i], NULL);
    }

    pthread_mutex_lock(&mutexsum);
    ActiveThreads -= 2;
    pthread_mutex_unlock(&mutexsum);

  }
}

void *threads_merge(void *thread_data){

  THREAD_DATA *tdata = thread_data;

  int tid = (*tdata).id; 
  int lo = (*tdata).lo;
  int cnt = (*tdata).cnt;
  int dir = (*tdata).dir;

  bitonicMerge(lo,cnt,dir);

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

      pthread_t threads[2];

      pthread_mutex_t mutexsum;

      pthread_mutex_init(&mutexsum, NULL);

      pthread_mutex_lock(&mutexsum);

      if (ActiveThreads < MaxThreads){

        ActiveThreads += 2;

        pthread_mutex_unlock(&mutexsum);
        
        THREAD_DATA t1;
        THREAD_DATA t2;
        
        t1.id = 0;
        t1.lo = lo;
        t1.cnt = cnt;
        t1.dir = ASCENDING;
  
        t2.id = 1;
        t2.lo = lo;
        t2.cnt = cnt;
        t2.dir = DESCENDING;

        pthread_create(&threads[0], NULL, threads_do, &t1);
        pthread_create(&threads[1], NULL, threads_do, &t2);

        parallel = 1;
      }
      
      if (parallel){

        int i;
        
        for (i=0; i<2; i++){
          pthread_join(threads[i], NULL);
        }


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

void *threads_do(void *thread_data){

  THREAD_DATA *tdata = thread_data;

  int tid = (*tdata).id;
  int lo = (*tdata).lo;
  int cnt = (*tdata).cnt;
  int dir = (*tdata).dir;

  int k = cnt/2;

  if(tid == 0){
    recBitonicSortParallel(lo,k,dir);
  }else{
    recBitonicSortParallel(lo+k,k,dir);
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

  int parallel = 0;

  pthread_t threads[2];

  pthread_mutex_t mutexsum;

  pthread_mutex_init(&mutexsum, NULL);

  pthread_mutex_lock(&mutexsum);

  if (ActiveThreads<MaxThreads){

    ActiveThreads += 2;

    pthread_mutex_unlock(&mutexsum);

    THREAD_DATA_IMP t1;
    THREAD_DATA_IMP t2;

    t1.id = 0;
    t1.start = start;
    t1.end = start+(end-start)/2;

    t2.id = 1;
    t2.start = start+(end-start)/2;
    t2.end = end;

    pthread_create(&threads[0], NULL, thread_imp, &t1);
    pthread_create(&threads[1], NULL, thread_imp, &t2);

    parallel = 1;
  }

  if (parallel){

    int l;

    for (l=0; l<2; l++){
      pthread_join(threads[l], NULL);
    }

    k = (end-start);

    for (j=k>>1; j>0; j=j>>1) {
      spawn_th(start,end-start,k,j);
      
    }

    pthread_mutex_lock(&mutexsum);
    ActiveThreads -= 2;
    pthread_mutex_unlock(&mutexsum);

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

void spawn_th(int start,int end,int k, int j){

  pthread_t threads[2];

  THREAD_DATA_IMP t1;
  THREAD_DATA_IMP t2;

  t1.start = start;
  t1.end = start+(end-start)/2;
  t1.k = k;
  t1.j = j;

  t2.start = start+(end-start)/2;
  t2.end = end;
  t2.k = k;
  t2.j = j;

  pthread_create(&threads[0], NULL, thread_imp_for, &t1);
  pthread_create(&threads[1], NULL, thread_imp_for, &t2);

  int l;

  for (l=0; l<2; l++){
    pthread_join(threads[l], NULL);
  }



}

void *thread_imp_for(void *t_data){

  THREAD_DATA_IMP *tdata = t_data;
  
  int start = (*tdata).start;
  int end = (*tdata).end;
  int k = (*tdata).k;
  int j = (*tdata).j;
  int i;

  for (i=start; i<end; i++) {
    int ij=i^j;
    if ((ij)>i) {
      if ((i&k)==0 && a[i] > a[ij]) 
        exchange1(i,ij);
      if ((i&k)!=0 && a[i] < a[ij])
        exchange1(i,ij);
    }
  }

  

}


void *thread_imp(void *t_data){

  THREAD_DATA_IMP *tdata = t_data;
  
  int tid = (*tdata).id;
  int start = (*tdata).start;
  int end = (*tdata).end;

  if (tid == 0){

    spawn(start,end);

  }else{
    qsort(a+start,end-start,sizeof(int),cmpfuncD);
  }

}

