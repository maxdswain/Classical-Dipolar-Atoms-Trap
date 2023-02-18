CC = gcc
CFLAGS = -Wall -g -Ofast -march=native
LDLIBS = -ltoml -lgsl -lgslcblas -lm

all: main reblock density

%: src/%.c src/shared.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm main reblock density *.txt *.out *.png
