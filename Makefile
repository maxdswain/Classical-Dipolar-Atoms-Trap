CC = gcc
CFLAGS = -Wall -g -O3
LDLIBS = -ltoml -lgsl -lgslcblas -lm

all: main reblock

main: src/main.c src/shared.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

reblock: src/reblock.c src/shared.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm main reblock *.txt *.out
