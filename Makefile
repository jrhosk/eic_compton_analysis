# Makefile to compile compton fitting code
#       Joshua Hoskins 
#         August 2017                                                                                                                             
#

ROOTLIBS   = $(shell root-config --libs ) -lSpectrum
ROOTGLIBS  = $(shell root-config --glibs)
LIB        = -L/usr/lib64/ -lboost_system -lboost_filesystem
INCLUDES   = -I$(shell root-config --incdir) -Iinclude/ -I/usr/include/
CC         = g++ ${INCLUDES}
SRC        = src
CFLAGS     = -O -Wall ${INCLUDES} ${LIB}

all: analysis

%.o: %.cc
	${CC} ${CFLAGS} -c -o $@ $< 
analysis : analysis.o ${SRC}/ComptonSimAnalysis.o ${SRC}/Generator.o ${SRC}/DataFile.o ${SRC}/vfTDCAnalysis.o
	${CC} ${INCLUDES} -o $@  ${CFLAGS} $^ ${ROOTLIBS} ${ROOTGLIBS} ${LIB}
clean:
	rm -f *.o *~ src/*.o src/*~ include/*~ config/*~
