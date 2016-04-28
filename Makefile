.PHONY: all clean

CCFLAGS= -I. -O3 -Wall
CPPCFLAGS=-I. -O3 -Wall

all: bin/parallel_motif

bin/parallel_motif: src/parallel_motif.cpp src/C-Thread-Pool/thpool.c
	mkdir -p bin
	gcc ${CCFLAGS} -c src/C-Thread-Pool/thpool.c -o src/thpool.o -lpthread
	mpic++ ${CPPCFLAGS} -c src/parallel_motif.cpp -o src/parallel_motif.o
	mpic++ src/thpool.o src/parallel_motif.o -o bin/parallel_motif
	rm src/*.o

clean:
	rm -r ./bin