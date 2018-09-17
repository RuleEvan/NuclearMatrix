CC=gcc
CFLAGS=-c -Wall -lm -ldl

all: wave

wave: main.o
	$(CC) main.o -o wave -lm -ldl -lgsl -lgslcblas

main.o: main.c
	$(CC) $(CFLAGS) main.c

clean:
	rm -rf *.o wave
