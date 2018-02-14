#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

int N,k,D;
double **dist; //Contains distances between each point and the k closest neighbors
double **ind;  //Contains the index's of the closest neighbors for each point 

struct timeval startwtime, endwtime;
double seq_time;

double search_time = 0.0;

MPI_Request reqs[2];   // required variable for non-blocking calls
MPI_Status stats[2];   // required variable for Waitall routine

MPI_Status status;

double **read_data_from_file(int rank,int numtasks);
double **alloc_2d_init(int rows, int cols);
double *init_L(int rank,int numtasks);
void knn_search(double **X,double **Y,int index,int numtasks,int rank,int first_time);
void insertionSort(double *arr, double *neighbors, int n);
void verify_results(double *Mnn,int *matches,int rank,int numtasks);
void makeCommunications(double **points,double **recv_points,int rank,int numtasks);
void updateTables(double **points_to_send, double **points_to_receive,int numtasks);
double **setLabelsOfClosestNeighbors(double **ind,int numtasks);
double *MostFrequentValuesOfClosestNeighbors(double **Lnn, int numtasks);


void main(int argc, char** argv) {

    if(argc != 4){
      printf("Invalid number of arguments! Aborting...");
      exit(0);
    }

    int rank;
    int numtasks;
    
    N = atoi(argv[1]);
    k = atoi(argv[2]);
    D = atoi(argv[3]);

    MPI_Init(&argc,&argv);

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double **points;
    double **recv_points;
    double **points_to_send;
    double **points_to_receive;

    double **Lnn;

    int first_time = 1;

    double *Mnn;
    int *matches;

    matches = (int *)malloc(N/numtasks*sizeof(double));

    dist = alloc_2d_init(N/numtasks,k);
    ind = alloc_2d_init(N/numtasks,k);

    Lnn = alloc_2d_init(N/numtasks,k);

    recv_points = alloc_2d_init(N/numtasks,D);
    points_to_send = alloc_2d_init(N/numtasks,D);
    points_to_receive = alloc_2d_init(N/numtasks,D);

    points = read_data_from_file(rank,numtasks);

    updateTables(points_to_send,points,numtasks);

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
      gettimeofday (&startwtime, NULL);

    makeCommunications(points_to_send,points_to_receive,rank,numtasks);

    knn_search(points,points,0,numtasks,rank,first_time);

    MPI_Waitall(2, reqs, stats);

    if(rank == 0){
      gettimeofday (&endwtime, NULL);

      seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
        + endwtime.tv_sec - startwtime.tv_sec);

      search_time += seq_time;
    }

    updateTables(recv_points,points_to_receive,numtasks);
    updateTables(points_to_send,recv_points,numtasks);

    first_time = 0;

    int i;
    int index = 0;
    int counter = 0;

    for(i=1; i<numtasks; i++){

      index = rank+i;

      if(index >= numtasks){
        index = counter;
        counter++;
      }
      
      if(rank == 0)
        gettimeofday (&startwtime, NULL);

      if(i<numtasks-1){
        makeCommunications(points_to_send,points_to_receive,rank,numtasks);
      }    
      
      knn_search(points,recv_points,index,numtasks,rank,first_time);

      if(i<numtasks-1){
        MPI_Waitall(2, reqs, stats);

        if(rank == 0){
          gettimeofday (&endwtime, NULL);

          seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
           + endwtime.tv_sec - startwtime.tv_sec);

          search_time += seq_time;
        }
        updateTables(recv_points,points_to_receive,numtasks);
        updateTables(points_to_send,recv_points,numtasks);
      }
    }
    
    if(rank == 0)
      printf("KNN SEARCH & COMMUNICATIONS CLOCK TIME = %f\n\n", search_time);

    Lnn = setLabelsOfClosestNeighbors(ind,numtasks);

    Mnn = MostFrequentValuesOfClosestNeighbors(Lnn,numtasks);

    verify_results(Mnn,matches,rank,numtasks);

    MPI_Finalize();
}

/*This is a function that copies one 2D array into another*/
void updateTables(double **points_to_send, double **points_to_receive,int numtasks){

  int i,j;

  for(i=0; i<N/numtasks; i++){
     for(j=0; j<D; j++){
        points_to_send[i][j] = points_to_receive[i][j];
      }
    }
}

/*This is a function which reads point coordinates from specified binary file*/
double **read_data_from_file(int rank,int numtasks){

  FILE *fp;
  double *buffer;
  double **ret_buf;
  int size,number_of_elements;
  int i,j,counter;
  int k;

  fp = fopen("X.bin","rb");
  
  fseek(fp, 0, SEEK_END);
  
  size = ftell(fp); //total size in bytes

  number_of_elements = size/8; //total number of elements

  buffer = (double *)malloc(N/numtasks*sizeof(double));

  ret_buf = alloc_2d_init(N/numtasks,D);

  fseek(fp,0,SEEK_SET);

  int index = rank*N/numtasks*8;

  for(i=0; i<D; i++){
    fseek(fp,index,SEEK_CUR);
    fread(buffer,sizeof(double),N/numtasks,fp);
    index = (numtasks-1)*N/numtasks*sizeof(double);

    for(k=0; k<N/numtasks; k++)
      ret_buf[k][i] = buffer[k];
  }

  fclose(fp);
  
  return ret_buf;
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

/*This is the function which executes KNN search on specified points*/
void knn_search(double **X,double **Y,int index,int numtasks,int rank,int first_time){

  int i,j,c;
  double distance;
  int counter;

  for(i=0; i<N/numtasks; i++){
    counter = 0;
    for(j=0; j<N/numtasks; j++){
      
      if(i!=j){
        
        distance = 0.0;
        for(c=0; c<D; c++){
          distance += (X[i][c]-Y[j][c])*(X[i][c]-Y[j][c]); 
        }

        if(counter<=k && first_time == 1){
          dist[i][counter] = distance;
          ind[i][counter] = j+rank*N/numtasks;
          counter++;
          if(counter == k)
            insertionSort(dist[i],ind[i],k);
        }else{
          if(distance<dist[i][k-1]){
            dist[i][k-1] = distance;
            if(first_time){
              ind[i][k-1] = j+rank*N/numtasks; 
            }else{
              ind[i][k-1] = j+index*N/numtasks;
            }       
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


/*This is the function that creates ring communications between processes*/
void makeCommunications(double **points,double **recv_points,int rank,int numtasks){
  
  int prev,next;

  if(rank == 0){
    prev = 1;
    next = numtasks-1;
  }else if(rank == numtasks-1){
    prev = 0;
    next = rank-1;
  }else{
    prev = rank+1;
    next = rank-1;
  }

  MPI_Irecv(&recv_points[0][0], N/numtasks*D, MPI_DOUBLE, prev, 1, MPI_COMM_WORLD, &reqs[0]);
  MPI_Isend(&points[0][0], N/numtasks*D, MPI_DOUBLE, next, 1, MPI_COMM_WORLD, &reqs[1]);

}
/*This is a function which creates a 2D array with the labels of the closest neighbors*/
double **setLabelsOfClosestNeighbors(double **ind,int numtasks){
  
  FILE *fp;
  int i,j;
  int key;
  double label;
  double **Lnn;
  int bytes_to_offset;

  Lnn = alloc_2d_init(N/numtasks,k);

  fp = fopen("L.bin","rb");

  for(i=0; i<N/numtasks; i++){
    for(j=0; j<k; j++){
      key = ind[i][j];
      bytes_to_offset = key*sizeof(double);
      fseek(fp,bytes_to_offset,SEEK_SET);
      fread(&label,sizeof(label),1,fp);
      Lnn[i][j] = label;
    }
  }

  fclose(fp);

  return Lnn;
}

/*This is a function which calculates the most frequent labels for each point*/
double *MostFrequentValuesOfClosestNeighbors(double **Lnn, int numtasks){
  
  double *Mnn;
  int i,j;
  double key;
  int count;
  int max;
  int pos_max;

  Mnn = (double *)malloc(N/numtasks*sizeof(double));

  for(i=0; i<N/numtasks; i++){
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
  return Mnn;
}

/*This is the function which calculates the percentage succeded*/
void verify_results(double *Mnn,int *matches,int rank,int numtasks){

  int i;
  int count;
  double *L;
  FILE *fp;

  L = init_L(rank,numtasks);
  
  for(i=0; i<N/numtasks; i++){
    if(L[i]-Mnn[i] == 0.0){
      matches[i] = 1;
    }else{
      matches[i] = 0;
    }
  }
  
  count = 0;
  for(i=0; i<N/numtasks; i++){
    if(matches[i] == 1){
      count++;
    }
  }

  printf("Percentage Succeded = %.6f\n", (double)count*numtasks/N*100);

}

/*This is a function that reads labels for each point int the process*/
double *init_L(int rank,int numtasks){

   FILE *fp;
   double *buffer;
   int size,number_of_elements;
   int i,j,counter;

   fp = fopen("L.bin","rb");

   fseek(fp, 0, SEEK_END);

   size = ftell(fp); //total size in bytes

   number_of_elements = size/8; //total number of elements

   buffer = (double *)malloc(number_of_elements/numtasks*sizeof(double));

   fseek(fp, rank*N/numtasks*8, SEEK_SET);

   fread(buffer,sizeof(buffer),number_of_elements/numtasks,fp);

   return buffer;
}
