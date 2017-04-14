#include "pgm.h"
#include "filters.h"
#include <math.h>
/*
    Filtrage linéaire résursif
    S'utilise sous la forme "filtreAdaptatifRecursif source resultat k"
 */

int main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
    int nl,nc;
    unsigned char ** im1=NULL;
    double** im4;

    if (ac < 6) {
        printf("Usage : %s entree sortie t r sigma\n", av[0]); exit(1);
    }
    /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
    im1=lectureimagepgm(av[1],&nl,&nc);
    if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }
    /* Calcul de son inverse video */
    double**im3=imuchar2double(im1,nl,nc);
    im4=alloue_image_double(nl,nc);

    //Filtrage
    NIMeansFilter(im3, im4, nl, nc, atoi(av[3]), atoi(av[4]), atof(av[5]));

    ecritureimagepgm(av[2],crop(imdouble2uchar(im4,nl,nc),0,0,nl,nc),nl,nc);
}
