CFLAGS=-c -g -std=gnu99
LDFLAGS=-lm

all: filtreBilateral

filtreGaussien: filtreBilateral.o pgm.o filtre.o 
	gcc $^ -o $@ $(LDFLAGS)

%.o:%.c
	gcc $(CFLAGS) $<

clean:
	rm -f *.o filtreBilateral
