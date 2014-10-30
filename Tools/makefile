
CC=gcc
CPP=gcc
CFLAGS=-O2 -march=x86-64 -fPIC -pipe
#CFLAGS=-fPIC -Wall
LDFLAGS=-lm



###
#
# commands to make 
#
###

cuda: allall fitzhugh_cuda leech_cuda
	rm *.o

all: allall
	rm *.o

###
#
# Link dynamic libraries
# (needs a directory ./lib/ which is created)
#
###

allall: tools fitzhugh thetax2 leech plant


tools: tools.o lib
	gcc -shared tools.o $(LDFLAGS) -o lib/_tools.so

thetax2: thetax2.o lib
	gcc -shared thetax2.o -o lib/_thetax2.so -lm

fitzhugh: fitzhugh.o lib
	gcc -shared fitzhugh.o -o lib/_fitzhugh.so -lm

plant: plant.o lib
	gcc -shared plant.o -o lib/_plant.so -lm

leech: leech.o lib
	gcc -shared leech.o -o lib/_leech.so -lm

leech_cuda: leech.o leech_cuda.o lib
	gcc -shared leech_cuda.o leech.o -o lib/_leech_cuda.so -L/usr/local/cuda/lib64 -lm -lcudart -lcuda

fitzhugh_cuda: fitzhugh.o fitzhugh_cuda.o lib
	gcc -shared fitzhugh_cuda.o fitzhugh.o -o lib/_fitzhugh_cuda.so -L/usr/local/cuda/lib64 -lm -lcudart -lcuda


###
#
# Compile programs
#
###


leech_cuda.o: leech.cu
	nvcc -c leech.cu -o leech_cuda.o -arch sm_20 -Xcompiler -fPIC

fitzhugh_cuda.o: fitzhugh.cu
	nvcc -c fitzhugh.cu -o fitzhugh_cuda.o -arch sm_20 -Xcompiler -fPIC

.c.o:
	$(CC) -c $< $(CFLAGS)

lib:
	mkdir lib


###


clean:
	rm -rf *.o *.pyc lib
