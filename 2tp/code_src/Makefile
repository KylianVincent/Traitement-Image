CFLAGS=-c -g -std=gnu99
LDFLAGS=-lm

all: filtreBilateral filtreAdaptatifRecursif filtreMedian filtreNLMeans estimationBruit testPsnr

filtreBilateral: filtreBilateral.o pgm.o filters.o
	gcc $^ -o $@ $(LDFLAGS)

filtreAdaptatifRecursif: filtreAdaptatifRecursif.o pgm.o filters.o
	gcc $^ -o $@ $(LDFLAGS)

filtreMedian: filtreMedian.o pgm.o filters.o
	gcc $^ -o $@ $(LDFLAGS)

filtreNLMeans: filtreNLMeans.o pgm.o filters.o
	gcc $^ -o $@ $(LDFLAGS)

estimationBruit: estimationBruit.o pgm.o filters.o
	gcc $^ -o $@ $(LDFLAGS)

testPsnr: testPsnr.o pgm.o
	gcc $^ -o $@ $(LDFLAGS)

%.o:%.c
	gcc $(CFLAGS) $<

clean:
	rm -f *.o filtreBilateral filtreAdaptatifRecursif filtreMedian filtreNLMeans estimationBruit testPsnr
