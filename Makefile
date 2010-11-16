# Makefile for GridrecIDL for Linux

SHELL = /bin/sh
CFLAGS = -O -DIDL -fPIC -Wall
INCLUDES=-I../
MAKE = make

# For Numerical Recipes
#OBJ = grid.o GridrecIDL.o grid_math.o pswf.o filters.o fft.o
#SRC = grid.c GridrecIDL.c grid_math.c pswf.c filters.c fft.c

# For FFTW
OBJ = grid.o GridrecIDL.o grid_math.o pswf.o filters.o fft_fftw.o
SRC = grid.c GridrecIDL.c grid_math.c pswf.c filters.c fft_fftw.c

GridrecIDL.so: $(OBJ)
	@echo "Creating GridrecIDL.so"
# For Numerical Recipes
#	$(LD) -G -o$@ $(OBJ) -lm;
# For FFTW
	$(LD) -G -o$@ $(OBJ) -lfftw3f -lm;

clean:
	rm *.o GridrecIDL.so 

