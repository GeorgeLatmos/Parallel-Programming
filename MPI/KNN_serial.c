#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

int N,k,D;
double **dist; //Contains distances between each point and the k closest neighbors
double **ind;  //Contains the index's of the closest neighbors for each point

struct timeval startwtime, endwtime;
double seq_time;

double **init();

double *init_L();
double **alloc_2d_init(int rows, int cols);
float **alloc_2d_init2(int rows, int cols);
void knn_search(double **X,double **Y);
void insertionSort(double *arr, double *neighbors, int n);
void set_Lnn(double **Lnn,double *L);
void set_Mnn(double *Mnn,double **Lnn);
void verify_results(double *Mnn, double *L,int *matches);


void main(int argc, char *argv[]) {

   if (argc != 4){
      printf("Invalid number of arguments! Aborting...\n");
      exit(0);
   }

   N = atoi(argv[1]);
   k = atoi(argv[2]);
   D = atoi(argv[3]);

   double **points; //Contains total points and number of coordinates
   double **Lnn;    //Contains the labels of the closest neighbors

   double *L;     //Defines a label for each point
   double Mnn[N];   //Find most frequent values in Lnn for each point
   int matches[N];
   
   L = (double *)malloc(N*sizeof(double));
   dist = alloc_2d_init(N,k);
   ind = alloc_2d_init(N,k);
   Lnn = alloc_2d_init(N,k);

   printf("......READ POINTS FROM FILE......\n");

   points = init();

   printf("......DONE......\n\n");

   L = init_L(); 

   printf("......KNN SEARCH......\n");

   gettimeofday (&startwtime, NULL);

   knn_search(points,points);

   gettimeofday (&endwtime, NULL);

   seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
          + endwtime.tv_sec - startwtime.tv_sec);

   printf("KNN SERIAL CLOCK TIME = %f\n\n", seq_time);

   printf("......MATCH LABELS......\n");

   gettimeofday (&startwtime, NULL);

   set_Lnn(Lnn,L);

   set_Mnn(Mnn,Lnn);

   gettimeofday (&endwtime, NULL);

   seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
          + endwtime.tv_sec - startwtime.tv_sec);

   printf("MATCHING LABELS CLOCK TIME = %f\n\n", seq_time);

   printf("......VERIFYING RESULTS......\n");
   
   gettimeofday (&startwtime, NULL);

   verify_results(Mnn,L,matches);

   gettimeofday (&endwtime, NULL);

   seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
          + endwtime.tv_sec - startwtime.tv_sec);

   printf("VERIFYING RESULTS CLOCK TIME = %f\n", seq_time);

}

double **alloc_2d_init(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    int i;
    for (i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}


double **init(){

   FILE *fp;
   double *buffer;
   double **ret_buf;
   int size,number_of_elements;
   int i,j,counter;

   fp = fopen("X.bin","rb");

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

double *init_L(){

   FILE *fp;
   double *buffer;
   int size,number_of_elements;
   int i,j,counter;

   fp = fopen("L.bin","rb");

   fseek(fp, 0, SEEK_END);

   size = ftell(fp); //total size in bytes

   number_of_elements = size/8; //total number of elements

   buffer = (double *)malloc(number_of_elements*sizeof(double));

   fseek(fp, 0, SEEK_SET);

   fread(buffer,sizeof(buffer),number_of_elements,fp);

   return buffer;
}

void knn_search(double **X,double **Y){

  int i,j,c;
  float distance;
  int counter;

  for(i=0; i<N; i++){
    counter = 0;
    for(j=0; j<N; j++){
      
      if(i!=j){
        
        distance = 0.0;
        for(c=0; c<D; c++){
          distance += (X[i][c]-Y[j][c])*(X[i][c]-Y[j][c]); 
        }

        if(counter<=k){
          dist[i][counter] = distance;
          ind[i][counter] = j;
          counter++;
          if(counter == k)
            insertionSort(dist[i],ind[i],k);
        }else{
          if(distance<dist[i][k-1]){
            dist[i][k-1] = distance;
            ind[i][k-1] = j;
            insertionSort(dist[i],ind[i],k);
          }
        }
      }
    }
  }
}

/*Function to sort an array using insertion sort*/
void insertionSort(double *arr, double *neighbors, int n)
{

   int i,j;
   double key,temp;

   for (i = 1; i < n; i++)
   {
       key = arr[i];
       j = i-1;
 
        /*Move elements of arr[0..i-1], that are
          greater than key, to one position ahead
          of their current position*/ 
       while (j >= 0 && arr[j] > key)
       {
           arr[j+1] = arr[j];
           temp = neighbors[j+1];
           neighbors[j+1] = neighbors[j];
           neighbors[j] = temp;
           j = j-1;
       }
       arr[j+1] = key;
   }
}

void set_Lnn(double **Lnn,double *L){
  
  int i,j;
  double key;

  for(i=0; i<N; i++){
    for(j=0; j<k; j++){
      key = ind[i][j];
      Lnn[i][j] = L[(int)key];
    }
  }

}

void set_Mnn(double *Mnn,double **Lnn){
  int i,j;
  double key;
  int count;
  int max;
  int pos_max;

  for(i=0; i<N; i++){
    max = -1;
    key = 0.0;
    while(key<(double)10){
      count = 0;
      for(j=0; j<k; j++){
        if(key == Lnn[i][j]){
          count++;
        }
      }
      if(count>max){
        max = count;
        pos_max = key;
      }
      key++;
    }
    Mnn[i] = (double)pos_max;
  }

}

void verify_results(double *Mnn, double *L,int *matches){

  int i;
  int count;
  
  for(i=0; i<N; i++){
    if(L[i]-Mnn[i] == 0.0){
      matches[i] = 1;
    }else{
      matches[i] = 0;
    }
  }
  
  count = 0;
  for(i=0; i<N; i++){
    if(matches[i] == 1){
      count++;
    }
  }

  printf("Percentage Succeded = %.6f\n", (double)count/N*100);

}
