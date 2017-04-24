CFLAGS=-c -g -std=gnu99
LDFLAGS=-lm

all: filtreBilateral filtreAdaptatifRecursif filtreMedian filtreNIMeans estimationBruit

filtreBilateral: filtreBilateral.o pgm.o filters.o
	gcc $^ -o $@ $(LDFLAGS)

filtreAdaptatifRecursif: filtreAdaptatifRecursif.o pgm.o filters.o
	gcc $^ -o $@ $(LDFLAGS)

filtreMedian: filtreMedian.o pgm.o filters.o
	gcc $^ -o $@ $(LDFLAGS)

filtreNIMeans: filtreNIMeans.o pgm.o filters.o
	gcc $^ -o $@ $(LDFLAGS)

estimationBruit: estimationBruit.o pgm.o filters.o
	gcc $^ -o $@ $(LDFLAGS)

%.o:%.c
	gcc $(CFLAGS) $<

clean:
	rm -f *.o filtreBilateral filtreAdaptatifRecursif filtreMedian filtreNIMeans estimationBruit
