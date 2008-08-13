# Makefile for fft_test for Linux

SHELL = /bin/sh
CFLAGS = -O2 -shared -DLINUX -L/opt/intel/mkl/lib/32
OBJ = fft_test.o ../fft.o
SRC = fft_test.c

fft_test.so: $(OBJ)
# If have MKL
#	$(CC) $(CFLAGS) -o $@ $(OBJ) -lguide -lmkl_def -lfftw3 -lm;
# If no MKL
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lfftw3f -lm;

clean:
	rm *.o fft_test.so 

