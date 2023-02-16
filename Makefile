CC = gcc
CFLAGS = -Wall -g -O3
LDLIBS = -ltoml -lgsl -lgslcblas -lm

all: main reblock density

%: src/%.c src/shared.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm main reblock *.txt *.out
