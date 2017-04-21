CC = gcc
CC_FLAGS =  -w -g -O2 -ftree-vectorize -fopt-info-vec-optimized -fopenmp
LFLAGS = -g -lhdf5 -lm -lfftw3 -lfftw3_threads -lgomp


all:

	$(CC) $(CC_FLAGS) -c main.c hdf5.c perf.c grid.c
	$(CC) -o gridder.out main.o hdf5.o perf.o grid.o $(LFLAGS)

test: 
	./gridder.out --theta=0.08 --lambda=300000 --image=image.out --grid=grid.out --wkern=../../data/kernels/SKA1_Low_wkern.h5 ../../data/vis/SKA1_Low_quick.h5
