#include "mpi.h"
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<complex.h>
#include<time.h>
#define MSIZE 10000
#define MAX_ITS 1000
#define C_REAL_MAX 1
#define C_REAL_MIN -2
#define C_IMG_MAX 1.5
#define C_IMG_MIN -1.5
#define C_RRANGE C_REAL_MAX - C_REAL_MIN
#define C_IRANGE C_IMG_MAX - C_IMG_MIN
#define TIMEREP 3 //overall time reported in 10^-(TIMEREP) seconds

int cal_pixel(double complex c)
{
 int count = 0;
 double temp, lengthsq;

 double complex z = 0;

 do
 {
  temp = creal(z) * creal(z) - cimag(z) * cimag(z) + creal(c);
  z = temp + (2*creal(z)*cimag(z) + cimag(c)) * I;
  lengthsq = creal(z) * creal(z) + cimag(z) * cimag(z);
  count++;
 }
 while((lengthsq < 4.0) && (count < MAX_ITS));

 if(count == MAX_ITS)
  return 1;
 else return 0;

}

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
 int nl, start, stop;
 int bsize, method = 0;
 int i, j;
 int sent, recd, redlight;
 int **mset, **mset2;
 int *data, *data2, *recbuf, *sendbuf;
 MPI_Status *deets;
 FILE *fp, *fp2;
 struct timespec tstart, tstop;

 MPI_Init(&argc, &argv);

 stat = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
 assert(stat == MPI_SUCCESS);

 stat = MPI_Comm_rank(MPI_COMM_WORLD, &mype);
 assert(stat == MPI_SUCCESS);


 if(argc > 2)
 {
  if(mype == 0)
   fprintf(stderr, "usage: %s [chunk size]\n", argv[0]);
  MPI_Finalize();
  exit(1);
 }

 /*check if the default method, or the master-slave method is being used
   if master-slave, set the block size to argv[1]*/
 if(argc == 2)
 {
  method = 1;
  bsize = atoi(argv[1]);
 }

 if(mype == 0)
 {
  mset2 = (int **)malloc(sizeof(int *)*MSIZE);
  data2 = (int *)malloc(sizeof(int)*MSIZE*MSIZE);

  for(i = 0; i < MSIZE; ++i)
  {
   mset2[i] = &data2[i*MSIZE];
  }
 }

 /*fixed distribution of work*/
 if(!method)
 {
  /*arrays used for displacement and sizes in allgatherv*/
  int displ[nprocs]; 
  int rcounts[nprocs];

  clock_gettime(CLOCK_MONOTONIC, &tstart);

  /*each process does nl lines of the final matrix */
  /*if MSIZE % nprocs != 0, add a line for each of the first MSIZE % nprocs processes*/
  nl = MSIZE/nprocs;
  start = nl*mype;
 
  for(i = 0; i < MSIZE % nprocs; ++i)
  {
   if(i == mype)
   {
    ++nl;
    start += i;
   }
  }

  for(i = MSIZE % nprocs; i < nprocs; ++i)
  {
   if(i == mype)
    start += (MSIZE % nprocs);
  }
  stop = start+nl;

  /*initializing buffers to be used for calculations, and a 
    2D array overlay for clarity*/
  mset = (int **)malloc(sizeof(int *)*nl);
  data = (int *)malloc(sizeof(int)*MSIZE*nl);

  for(i = 0; i < nl; ++i)
  {
   mset[i] = &data[i*MSIZE];
  }

  
  for(i = start; i < stop; ++i)
  {
   for(j = 0; j < MSIZE; ++j)
    mset[i-start][j] = cal_pixel(((double)C_REAL_MIN + (((double)C_RRANGE)/MSIZE)*j) +  ((double)C_IMG_MAX - (((double)C_IRANGE)/MSIZE)*i)*I);
  }

  /*Accounting for non-uniform sizes and displacements in the event that MSIZE % nprocs != 0*/
  for(i = 0; i < nprocs; ++i)
  {
   rcounts[i] = (MSIZE/nprocs) * MSIZE;
   displ[i] = i*(MSIZE/nprocs)*MSIZE;
  }
  
  for(i = 0; i < MSIZE % nprocs; ++i)
  {
   rcounts[i] += MSIZE;
   displ[i] += i*MSIZE;
  }

  for(i = MSIZE % nprocs; i < nprocs; ++i)
   displ[i] += MSIZE * (MSIZE % nprocs);
 
  MPI_Gatherv(data, nl*MSIZE, MPI_INT, data2, rcounts, displ, MPI_INT, 0, MPI_COMM_WORLD);

  free(data);
  free(mset);
  clock_gettime(CLOCK_MONOTONIC, &tstop);
 }

 /*master/slave process code*/
 else
 {
  if(nprocs < 2)
  { 
   if(mype == 0)
    fprintf(stderr, "Master/slave process requires at least 2 processors\n");
   MPI_Finalize();
   exit(1);
  }
 
  clock_gettime(CLOCK_MONOTONIC, &tstart);
  
  deets = (MPI_Status *)malloc(sizeof(MPI_Status));

  if(mype == 0)
  {
   recbuf = (int *)malloc(sizeof(int)*bsize);
   
   sent = bsize*(nprocs - 1); 	//Number of points already sent out to be calculated
   recd = 0;			//Number of points received back after calculation
   stop =  (MSIZE*MSIZE)/bsize + (((MSIZE*MSIZE) % bsize == 0) ? 0 : 1); //Number of total receipts expected
   redlight = -1;		//Message to send out to tell processes to stop waiting for more work

   while(recd < stop)
   {
    MPI_Recv(recbuf, bsize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, deets);

    start = deets->MPI_TAG;	//tag is used as an indication of where in the matrix of points this block should start
    for(i = start; (i < (start + bsize)) && (i < MSIZE*MSIZE); ++i)
     mset2[i / MSIZE][i % MSIZE] = recbuf[i - start];
    ++recd;
    /*if not all of the points have been sent out for calculation, send the next chunk to the 
     node just received from*/
    if(sent < MSIZE*MSIZE)
    {
     MPI_Send(&sent, 1, MPI_INT, deets->MPI_SOURCE, 0, MPI_COMM_WORLD);
     sent += bsize;
    }
    /*otherwise, send that node a redlight*/
    else
    { 
     MPI_Send(&redlight, 1, MPI_INT, deets->MPI_SOURCE, 0, MPI_COMM_WORLD);
    }
   }
   free(recbuf);
  }
 
  else
  {
   /*start by having each node do 1 chunk of calculations in order*/
   start = (mype-1) * bsize;
   sendbuf = (int *)malloc(sizeof(int)*bsize);
  
   /*wait for master node to send the number of the starting node for the next chunk of calculations
    if a redlight is received instead, stop calculating*/
   /*if nprocs * bsize > the total number of calculations, some nodes may initially try to calculate 
     points that are outside the scope of the matrix. Hence the second parameter for the loop*/
   while((start != -1) && (start < (MSIZE*MSIZE)))
   {
    for(i = start; (i < (start + bsize)) && (i < (MSIZE*MSIZE)); ++i)
     sendbuf[i-start] = cal_pixel(((double)C_REAL_MIN + (((double)C_RRANGE)/MSIZE)*(i%MSIZE)) + 
 	((double)C_IMG_MAX - (((double)C_IRANGE)/MSIZE)*(i/MSIZE))*I);
    MPI_Send(sendbuf, bsize, MPI_INT, 0, start, MPI_COMM_WORLD);
    MPI_Recv(&start, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, deets);
   }
   
   free(sendbuf);
  }

 free(deets);


 clock_gettime(CLOCK_MONOTONIC, &tstop);
 }

 /*data print out*/
 if(mype == 0)
 {
  fp = fopen("dots.dat", "w");
  for(i = 0; i < MSIZE; ++i)
  {
   for(j = 0; j < MSIZE; ++j)
    fprintf(fp, "%d ", mset2[i][j]);
   fprintf(fp, "\n");
  }
  fclose(fp);
 
  fp2 = fopen(((method)? "ms-times" : "def-times"), "a");
  fprintf(fp2, "%d\t%.2lf\n", nprocs, calc_time(tstart, tstop, TIMEREP));
  fclose(fp2);
 } 


 /*cleanup*/

 if(mype == 0)
 {
  free(mset2[0]);
  free(mset2);
 }

 MPI_Finalize();
 exit(0);
}
