Asad Ali
MPCS 51087
Homework 2

Exercise 1
**********
Files:
 ex1.c : source code
 makefile: make file
 m-def.sbatch: sbatch for running default (static distribution) algorithm on 2, 4, 8, 16, and 32 procs
 m-ms.sbatch: sbatch for running master-slave (dynamic distribution) algorithm on 2, 4, 8, 16, and 32 procs
 
Compilation:
 $module load openmp
 $make [CFLAGS] [OUTPUT]
  Default:
   CFLAGS=-Wall -g
   OUTPUT=mandelbrot


 $make clean
 removes the output file, for what that's worth.

Running:
 $mpiexec -n nprocs ./mandelbrot [chunk size]

 Without a specified chunk size, mandelbrot runs in default (statically allocated per-process workload) mode. With a chunk size, it runs with process 0 as the master, dynamically distributing chunk size blocks of work to the other processes.

Output:
 If no chunk size is specified, the number of processes and the time it took to complete calculation is appended to a file called "def-times". Otherwise, it is appended to a file called "ms-times". In both cases, the set itself is output to a file called "dots.dat", which is a matrix of bits that act as a pixel map.

Results:
 These are found in the "results" folder. "plot.pg" is a script to plot
def-times and ms-times on a logarithmic scale using gnuplot (output to
"data.png"). "man-make.pg" is a script to plot out the actual Mandelbrot set
that was calculated, and outputs to "man-map.png". From what I could see,
there are significant order-of-magnitude speedups as the number of processors
doubles, but those start to taper off by 32 processors for both types of
algorithm. I would assume that, at that point, the message passing starts to
become a more significant portion of processing time. Also, as expected, the
master-slave algorithm sped up faster, although it started off slower. That
makes sense, since, for two processes, the default algorithm has both actually
doing work, while the m-s algorithm relies on just the one.

Exercise 2
**********
Files:
 ex2.c: source code
 makefile: make file
 btes.sbatch: sbatch to run tests on 2 nodes
 
Compilation:
 $module load openmpi
 $make [CFLAGS] [OUTPUT]
  Default:
   CFLAGS=-Wall -g
   OUTPUT=btes

 $make clean
  removes btes

Running:
 $mpiexec -n 2 ./btes

 I suppose you could run it with more than 2 processes, but that just seems silly. As such, my program does not stand for it.

Output:
 Latency and bandwidth measurements are output to a file called "data". Latency is put as a comment at the top. Below it, two columns show the size of
the message passed (in bytes), and the bandwidth measured (in bytes/sec). 

Results:
 Again, these are found in the "results" folder. I made two plotting scripts -
one linear, and one logarithmic. Both show a sharp increase in bandwidth that
plateaus around 1GB/s at about the 256KB message size range, which leads me to believe that that's about the saturation point for bandwidth. The linear plot shows the plateau more starkly while the logarithmic plot lets you see all the data more clearly. 
 The latency measurement I got was consistently around 0.0009ms.
