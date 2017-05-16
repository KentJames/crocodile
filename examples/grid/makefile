CC = gcc
CC_FLAGS =  -w -g -O3 -ftree-vectorize -fopt-info-vec-optimized -ffast-math -fopenmp
CFLAGS = 
LFLAGS = -g -lhdf5 -lm -lfftw3 -lfftw3_threads -lgomp


all: grid

grid: grid.c main.c hdf5.c perf.c


	$(CC) $(CFLAGS) $(CC_FLAGS) -c main.c hdf5.c perf.c
	$(CC) $(CFLAGS) $(CC_FLAGS) -c grid.c 
	$(CC) -o gridder.out main.o hdf5.o perf.o grid.o $(LFLAGS)

test: 
	./gridder.out --threads=8 --theta=0.08 --lambda=300000 --wkern=../../data/kernels/SKA1_Low_wkern.h5 ../../data/vis/SKA1_Low_vis.h5
