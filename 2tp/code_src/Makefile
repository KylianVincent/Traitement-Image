CFLAGS=-c -g -std=gnu99
LDFLAGS=-lm

all: essai1 essai2 filtreGaussien detectionContourGrad

test: testPsnr filtreGaussien
	./test.sh

essai1: essai1.o pgm.o
	gcc $^ -o $@ $(LDFLAGS)

essai2: essai2.o pgm.o fft.o
	gcc $^ -o $@ $(LDFLAGS)

filtreGaussien: filtreGaussien.o pgm.o fft.o img.o
	gcc $^ -o $@ $(LDFLAGS)

testPsnr: testPsnr.o pgm.o
	gcc $^ -o $@ $(LDFLAGS)

detectionContourGrad: detectionContourGrad.o pgm.o img.o
	gcc $^ -o $@ $(LDFLAGS)

detectionContourLapl: detectionContourLapl.o pgm.o fft.o img.o
	gcc $^ -o $@ $(LDFLAGS)

diffTimefiltreGaussien: diffTimefiltreGaussien.o pgm.o fft.o img.o
	gcc $^ -o $@ $(LDFLAGS)

%.o:%.c
	gcc $(CFLAGS) $<

clean:
	rm -f *.o essai1 essai2 filtreGaussien detectionContourGrad testPsnr detectionContourLapl diffTimefiltreGaussien
