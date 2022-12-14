CC = gcc
CFLAGS = -Wall -g -O3
LDLIBS = -ltoml -lgsl -lgslcblas -lm

all: main

main: src/main.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm main *.txt
