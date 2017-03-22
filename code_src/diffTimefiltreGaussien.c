#include "pgm.h"
#include "img.h"
#include <math.h>
#include <time.h>
/*
    Debruitage par la methode de gauss
    S'utilise sous la forme  "exemple tangram.pgm res.pgm sigma (W)"
 */

int main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
    int nl,nc, oldnl,oldnc;
    int w;
    unsigned char ** im1=NULL;
    double** im4,** im41,** im5, ** im6, ** im7, **im8;

    if (ac < 2) {
        printf("Usage : %s entree",av[0]); exit(1);
    }
    /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
    im1=lectureimagepgm(av[1],&nl,&nc);
    if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }
    /* Calcul de son inverse video */
    double**im3=imuchar2double(im1,nl,nc);
    oldnl=nl; oldnc=nc;
    im4=alloue_image_double(nl,nc);

    double**im31=imuchar2double(im1,nl,nc);
    oldnl=nl; oldnc=nc;
    im41=alloue_image_double(nl,nc);

    for (int sigma = 0; sigma<100; sigma++) {
        clock_t debutFreq=clock();
        //FREQUENTIEL
        im7=padimdforfft(im3,&nl,&nc);
        im5=alloue_image_double(nl,nc); im6=alloue_image_double(nl,nc);
        fft(im7,im4,im5,im6,nl,nc);
        im8=alloue_image_double(nl,nc);
        /* Application du filtre gaussien */
        fftshift(im5, im6, im7, im8, nl, nc);
        filtrageFrequentiel(im7, im8, sigma/10.0, nl, nc);
        fftshift(im7, im8, im5, im6, nl, nc);
        ifft(im5,im6,im7,im8,nl,nc);
        clock_t finFreq=clock();

        clock_t debutSpace=clock();
        //SPATIAL
        w = ceil(3.0*sigma/10.0);
        filtrageSpatial(im31, im41, sigma/10.0, nl, nc, w);
        clock_t finSpace=clock();

        printf("%f\n", 1000*((double)(debutFreq-finFreq)-(debutSpace-finSpace))/CLOCKS_PER_SEC);
    }
}
