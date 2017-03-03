CFLAGS=-c -g -std=gnu99
LDFLAGS=-lm

all: essai1 essai2 filtreGaussien detectionContourGrad

test: testPsnr filtreGaussien
	./test.sh

essai1: essai1.o pgm.o
	gcc $(LDFLAGS) $^ -o $@

essai2: essai2.o pgm.o fft.o
	gcc $(LDFLAGS) $^ -o $@

filtreGaussien: filtreGaussien.o pgm.o fft.o
	gcc $(LDFLAGS) $^ -o $@

testPsnr: testPsnr.o pgm.o
	gcc $(LDFLAGS) $^ -o $@

detectionContourGrad: detectionContourGrad.o pgm.o
	gcc $(LDFLAGS) $^ -o $@

detectionContourLapl: detectionContourLapl.o pgm.o fft.o
	gcc $(LDFLAGS) $^ -o $@

%.o:%.c
	gcc $(CFLAGS) $<

clean:
	rm -f *.o essai1 essai2 filtreGaussien detectionContourGrad testPsnr detectionContourLapl
