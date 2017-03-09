#include "pgm.h"
#include "img.h"
#include <math.h>
/*
    Debruitage par la methode de gauss
    S'utilise sous la forme  "exemple tangram.pgm res.pgm sigma (W)"
 */

int main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
    int nl,nc, oldnl,oldnc;
    unsigned char ** im1=NULL;
    double** im4,** im5, ** im6, ** im7, **im8;

    if (ac < 4) {
        printf("Usage : %s entree sortie sigma (W)\n------ Si W est renseigné, un filtrage spatial est effectué, sinon on effectue un filtrage fréquentiel ------\n",av[0]); exit(1);
    }
    /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
    im1=lectureimagepgm(av[1],&nl,&nc);
    const double SIGMA = atof(av[3]);
    if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }
    /* Calcul de son inverse video */
    double**im3=imuchar2double(im1,nl,nc);
    oldnl=nl; oldnc=nc;
    im4=alloue_image_double(nl,nc);

    if (ac == 4) {
        // Utilisation du filtrage gaussien en fréquentiel
        /*  la fft demande des puissances de 2. On padde avec des 0, mais les dimensions nl et nc changent */
        im7=padimdforfft(im3,&nl,&nc);
        /*
          On peut aussi directement utiliser
          im7=padimucforfft(im1,&nl,&nc);
          sans convertir im1 en image de réels
        */
        /* Creation des images pour les parties reelles et imagianires des fft  */
        im5=alloue_image_double(nl,nc); im6=alloue_image_double(nl,nc);
        /* Calcul de la fft de im7,im4 */
        fft(im7,im4,im5,im6,nl,nc);
        /* Creation des images pour les parties reelles et imagianires des fft inverses */

        im8=alloue_image_double(nl,nc);
        /* Application du filtre gaussien */
        fftshift(im5, im6, im7, im8, nl, nc);
        filtrageFrequentiel(im7, im8, SIGMA, nl, nc);
        fftshift(im7, im8, im5, im6, nl, nc);

        /* Calcul de la fft inverse de im5,im6
           Le resultat doit etre stocke dans im4 */
        ifft(im5,im6,im4,im7,nl,nc);

    } else {
        // Utilisation du filtrage spatial
        const int W = atoi(av[4]);
        filtrageSpatial(im3, im4, SIGMA, nl, nc, W);
    }

    ecritureimagepgm(av[2],crop(imdouble2uchar(im4,nl,nc),0,0,oldnl,oldnc),oldnl,oldnc);
}
