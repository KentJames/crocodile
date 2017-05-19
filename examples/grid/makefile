CC = gcc
CC_FLAGS =  -w -g -O3 -ftree-vectorize -ffast-math -fopenmp
CFLAGS = 
LFLAGS = -g -lhdf5 -lm -lfftw3 -lfftw3_threads -lgomp


all: gridder.out 


clean: 
		rm *.o
		rm gridder.out

gridder.out: main.o hdf5.o perf.o grid.o
		$(CC) main.o hdf5.o perf.o grid.o -o gridder.out $(LFLAGS)


main.o: main.c
		$(CC) $(CFLAGS) $(CC_FLAGS) -c main.c

hdf5.o: hdf5.c 
		$(CC) $(CFLAGS) $(CC_FLAGS) -c hdf5.c

perf.o: perf.c 
		$(CC) $(CFLAGS) $(CC_FLAGS) -c perf.c

grid.o: grid.c
		$(CC) $(CFLAGS) $(CC_FLAGS) -c grid.c



test: 
	./gridder.out --threads=1 --theta=0.08 --lambda=300000 --image=image.out --wkern=../../data/kernels/SKA1_Low_wkern.h5 ../../data/vis/SKA1_Low_quick.h5

.PHONY: all clean test
