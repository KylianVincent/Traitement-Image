CFLAGS=-c -g -std=gnu99
LDFLAGS=-lm

all: filtreBilateral

filtreBilateral: filtreBilateral.o pgm.o filters.o 
	gcc $^ -o $@ $(LDFLAGS)

%.o:%.c
	gcc $(CFLAGS) $<

clean:
	rm -f *.o filtreBilateral
