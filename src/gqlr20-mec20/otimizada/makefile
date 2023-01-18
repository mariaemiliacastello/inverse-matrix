CFLAGS = -Wall -DLIKWID_PERFMON -I${LIKWID_INCLUDE}
LDFLAGS = -lm -L${LIKWID_LIB} -llikwid
DEBUG_FLAGS = -g
COMPILER = -std=gnu99

all: main

main: main.o sislin.o linhaComando.o utils.o
	gcc -o invmat main.o linhaComando.o sislin.o utils.o $(CFLAGS) $(LDFLAGS) $(COMPILER) $(DEBUG_FLAGS)

main.o: main.c
	gcc -c main.c -o main.o $(CFLAGS)

sislin.o: sislin.c utils.c sislin.h
	gcc -c sislin.c utils.c $(CFLAGS)

linhaComando.o: linhaComando.c linhaComando.h
	gcc -c linhaComando.c $(CFLAGS)	

utils.o: utils.c utils.h
	gcc -c utils.c $(CFLAGS)

clean:
	-rm -f *~ *.o invmat
