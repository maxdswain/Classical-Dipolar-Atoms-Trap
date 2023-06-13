CC = gcc
CFLAGS = -Wall -g -Ofast -march=native
LDLIBS = -ltoml -lgsl -lgslcblas -lm
BINS = main reblock density

all: $(BINS)

%: src/%.c src/shared.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm $(BINS) *.out *.png

.PHONY: all clean
