CC = gcc
CC_FLAGS = -Wall -g -O3 -ftree-vectorize -ffast-math -fopenmp
CFLAGS = 
LFLAGS = -g -lhdf5 -lm -lfftw3 -lfftw3_threads -lgomp


all: gridder.out 


clean: 
		rm *.o
		rm *.out

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

wtowers_test: 
	./gridder.out --theta=0.08 --lambda=300000 --wkern=../../data/kernels/SKA1_Low_wkern.h5 --subgrid=100 --margin=20 --winc=10 ../../data/vis/SKA1_Low_vis.h5

wtowers_test_vla: 
	./gridder.out --theta=0.1 --lambda=18000 --image=image.out --wkern=../../data/kernels/vlaa_wkern-150.h5 --subgrid=200 --margin=20 --winc=100 ../../data/vis/vlaa_theta0.1.h5

test: 
	./gridder.out --threads=12 --theta=0.08 --lambda=300000 --image=image.out --wkern=../../data/kernels/SKA1_Low_wkern.h5 ../../data/vis/SKA1_Low_vis.h5

test_vla:
	./gridder.out --theta=0.1 --lambda=18000 --image=image.out --wkern=../../data/kernels/vlaa_wkern-150.h5 ../../data/vis/vlaa_theta0.1.h5

.PHONY: all clean test
