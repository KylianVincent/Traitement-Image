#include "pgm.h"
#include "filters.h"
#include <math.h>
/*
    Estimation du bruit d'une image
    S'utilise sous la forme "estimationBruit image taille_blocs pourcentile"
 */

int main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
    int nl,nc;
    unsigned char ** im1=NULL;
    double** im2;

    if (ac < 4) {
        printf("Usage : %s image taille_blocs pourcentile\n", av[0]); exit(1);
    }
    /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
    im1=lectureimagepgm(av[1],&nl,&nc);
    if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }
    im2=imuchar2double(im1,nl,nc);

    //Estimation Bruit
    double varNoise;
    varNoise = noiseEstimation(im2, nl, nc, atoi(av[2]), atof(av[3]));

    printf("Variance du bruit = %f\n", varNoise);
}
