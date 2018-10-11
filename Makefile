CC=gcc
CFLAGS=-c -Wall -lm -ldl

all: nuclear

nuclear: main.o angular.o romberg.o phase_space.o potential.o brody.o av18.o matrix_element.o lanczos.o slater.o wave_function.o density.o master_eq.o 
	$(CC) main.o angular.o romberg.o phase_space.o potential.o brody.o av18.o matrix_element.o lanczos.o slater.o wave_function.o density.o master_eq.o -o nuclear -lm -ldl -lgsl -lgslcblas

main.o: main.c
	$(CC) $(CFLAGS) main.c

slater.o: slater.c
	$(CC) $(CFLAGS) slater.c

lanczos.o: lanczos.c
	$(CC) $(CFLAGS) lanczos.c

matrix_element.o: matrix_element.c
	$(CC) $(CFLAGS) matrix_element.c

potential.o: potential.c
	$(CC) $(CFLAGS) potential.c

av18.o: av18.c
	$(CC) $(CFLAGS) av18.c

brody.o: brody.c
	$(CC) $(CFLAGS) brody.c

phase_space.o: phase_space.c
	$(CC) $(CFLAGS) phase_space.c

romberg.o: romberg.c
	$(CC) $(CFLAGS) romberg.c

angular.o: angular.c
	$(CC) $(CFLAGS) angular.c

wave_function.o: wave_function.c
	$(CC) $(CFLAGS) wave_function.c

master_eq.o: master_eq.c
	$(CC) $(CFLAGS) master_eq.c

density.o: density.c
	$(CC) $(CFLAGS) density.c

clean:
	rm -rf *.o nuclear
