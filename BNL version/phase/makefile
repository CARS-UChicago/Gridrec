CFLAGS = -g  
CC = cc
#NETCDF = ../netcdf
HOME = /usr/people/jones
BIN = $(HOME)/bin
LIBS = -L$(HOME)/lib
INCLUDES = -I$(HOME)/include
LIBO = mema.o byteswap.o readPMIS.o readPXL.o shellsort.o fexist.o

all: PMIStoCDF dumpPMIS maxcPMIS pmislib.a pxllib.a phase1 phase2 nphase1 gvolume volume local $(LIBO)

PMIStoCDF: PMIStoCDF.c PMIStoCDF.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BIN)/$@ $@.c $(LIBS) -lnetcdf -lm;

dumpPMIS: dumpPMIS.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BIN)/$@ $@.c $(LIBS) -lnetcdf -lm

maxcPMIS: maxcPMIS.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BIN)/$@ $@.c $(LIBS) -lnetcdf -lm

mema.o: mema.c

readPMIS.o: readPMIS.c PMIStoCDF.h

byteswap.o: byteswap.c

shellsort.o: shellsort.c

fexist.o: fexist.c

pmislib.a: $(LIBO)
	ar rs pmislib.a $(LIBO)

phase1: phase1.c pmislib.a
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BIN)/$@ $@.c pmislib.a -lm

phase2: phase2.c 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BIN)/$@ $@.c pmislib.a -lm ${LIBS} -lnetcdf

pxllib.a: $(LIBO)
	ar rs pxllib.a $(LIBO)

nphase1: nphase1.c pxllib.a
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BIN)/$@ $@.c pxllib.a -lm

volume: volume.c 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BIN)/$@ $@.c $(LIBS) -lnetcdf -lm

gvolume: gvolume.c 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BIN)/$@ $@.c $(LIBS) -lnetcdf -lm

local: local.c 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BIN)/$@ $@.c $(LIBS) -lnetcdf -lm

clean:
	rm *.o
	rm *.a
