#include "pgm.h"
#include "img.h"
#include <math.h>

/*
    Detection de contours par filtre passe bande
    Utilisation de la TF du filtre LoG
 */

int main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
    int nl,nc, oldnl,oldnc;
    unsigned char ** im1=NULL;
    double** im4,** im5, ** im6, ** im7, **im8, **im9;

    if (ac < 4) {printf("Usage : %s entree sortie sigma\n",av[0]); exit(1); }
    im1=lectureimagepgm(av[1],&nl,&nc);
    const double SIGMA = atof(av[3]);

    if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }

    double**im3=imuchar2double(im1,nl,nc);
    oldnl=nl; oldnc=nc;
    im4=alloue_image_double(nl,nc);
    im4=padimdforfft(im3,&nl,&nc);
    /* Creation des images pour les parties reelles et imagianires des fft  */
    im5=alloue_image_double(nl,nc);
    im6=alloue_image_double(nl,nc);
    im7=alloue_image_double(nl,nc);
    /* Calcul de la fft de im4,im5 */
    fft(im4,im5,im6,im7,nl,nc);

    /* Creation des images pour les parties reelles et imagianires des fft inverses */
    im8=alloue_image_double(nl,nc); im9=alloue_image_double(nl,nc);
    /* Application du filtre gaussien */
    fftshift(im6, im7, im8, im9, nl, nc);
    filtrePasseBandeFrequentiel(im8, im9, nl, nc, SIGMA);
    //filtrePasseBandeIdeal(im8, im9, im6, im7, nl, nc);
    fftshift(im8, im9, im6, im7, nl, nc);

    /* Calcul de la fft inverse de im8,im9 */
    ifft(im6,im7,im8,im9,nl,nc);

    passageParZero(im8, im4, nl, nc);

    ecritureimagepgm(av[2],crop(imdouble2uchar(im4,nl,nc),0,0,oldnl,oldnc),oldnl,oldnc);
}
