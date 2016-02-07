#include "mpi.h"
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>
#define BASE_TEST 1000 //starting size for bandwidth test, in bytes
/*10 tests are run, with each one doubling the size of the last one*/
#define NUM_TESTS 13 //number of data points to plot for bandwidth
#define B_TEST 100 //Number of times to pass back and forth the bufer for each size in the bandwidth test
#define LAT_TEST 1000000 //number of times to pass back and forth the ping bufer for the latency test


/*takes two timespecs and returns a double representing the time difference
 *  in 10^-(power) seconds*/
double calc_time(struct timespec start, struct timespec stop, int power)
{
 int i;
 double seconds;
 double nseconds;
 
 seconds = stop.tv_sec - start.tv_sec;
 nseconds = stop.tv_nsec - start.tv_nsec;
 if(nseconds < 0)
 {
  nseconds += 1000000000;
  --seconds;
 }
 seconds = seconds + (nseconds / 1000000000);

 for(i = 0; i < power; ++i)
  seconds *= 10;

 return seconds;
}


int main(int argc, char **argv)
{
 int nprocs;
 int mype;
 int stat;
 int i, j, size = BASE_TEST;
 struct timespec start, stop;
 int **buf;
 int pingbuf[32];
 double latency;
 MPI_Status *deets;
 FILE *fp;

 MPI_Init(&argc, &argv);

 stat = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
 assert(stat == MPI_SUCCESS);

 stat = MPI_Comm_rank(MPI_COMM_WORLD, &mype);
 assert(stat == MPI_SUCCESS);

 if(nprocs != 2)
 {
  if(mype == 0)
   fprintf(stderr, "%s requires 2 processes\n", argv[0]);
  MPI_Finalize();
  exit(1);
 }

 for(i = 0; i < 32; ++i)
  pingbuf[i] = i;

 /*A series of buffers for the data points to be measured*/
 /*each buffer is double the size of the one before it*/
 buf = (int **)malloc(sizeof(int *) * NUM_TESTS);
 for(i = 0; i < NUM_TESTS; ++i)
 {
  buf[i] = (int *)malloc(sizeof(int)*size);
  for(j = 0; j < size; ++j)
   buf[i][j] = j;
  size *= 2;
 }

 deets = (MPI_Status *)malloc(sizeof(MPI_Status));

 /*latency test*/

 /*This consists of sending the "ping" buffer back and forth LAT_TEST times*/
 if(mype == 0)
 {
  fp = fopen("data", "w");
  clock_gettime(CLOCK_MONOTONIC, &start);
  for(i = 0; i < LAT_TEST; ++i)
  {
   MPI_Send(pingbuf, 32, MPI_INT, 1, 0, MPI_COMM_WORLD);
   MPI_Recv(pingbuf, 32, MPI_INT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, deets);
  }
  clock_gettime(CLOCK_MONOTONIC, &stop);

  /*The latency is measured as the overall time divided by the number of trips made*/
  latency = calc_time(start, stop, 3)/(2*LAT_TEST);
  fprintf(fp, "# Latency: %lfms\n", latency);
 }

 else if(mype == 1)
 {

  for(i = 0; i < LAT_TEST; ++i)
  {
   MPI_Recv(pingbuf, 32, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, deets);
   MPI_Send(pingbuf, 32, MPI_INT, 0, 0, MPI_COMM_WORLD);
  }

 }

 /*bandwidth tests*/
 /*For each of the NUM_TEST buffers, send the buffer back and forth between procs
   0 and 1 B_TEST times, then divide the total time by the number of trips made*/
 if(mype == 0)
 {
  size = BASE_TEST;
  for(i = 0; i < NUM_TESTS; ++i)
  {
   clock_gettime(CLOCK_MONOTONIC, &start);
   for(j = 0; j < B_TEST; ++j)
   {
    MPI_Send(buf[i], size, MPI_INT, 1, 0, MPI_COMM_WORLD);
    MPI_Recv(buf[i], size, MPI_INT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, deets);
   }
   clock_gettime(CLOCK_MONOTONIC, &stop);

   fprintf(fp, "%d\t%lf\n", size, 
	((double)size)/((calc_time(start, stop, 0) - (latency/1000))/(2*B_TEST)));
   size*=2;
  }
 }

 else if(mype == 1)
 {

  size = BASE_TEST;
  for(i = 0; i < NUM_TESTS; ++i)
  {
   for(j = 0; j < B_TEST; ++j)
   {
    MPI_Recv(buf[i], size, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, deets);
    MPI_Send(buf[i], size, MPI_INT, 0, 0, MPI_COMM_WORLD);
   }
   size *= 2;
  }

 }


 if(mype == 0)
  fclose(fp);

 MPI_Finalize();
 exit(0);
}
