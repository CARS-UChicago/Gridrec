# Makefile for grid and center 12:00 PM 10/17/97

SHELL = /bin/sh/
HOME=/usr/people/jones/
BIN=$(HOME)bin/
CFLAGS = -O2 
INCLUDES=-I$(HOME)include
MAKE = make
OBJ = grid.o gridmain.o grid_io.o pswf.o filters.o fft.c
SRC = grid.c gridmain.c grid_io.c pswf.c filters.c fft.c

gridrec: $(OBJ)
	@echo "Creating gridrec"
	$(CC) $(CFLAGS) -o $(BIN)$@ $(OBJ) -L$(HOME)lib -lnetcdf -lm;

clean:
	rm *.o gridrec

