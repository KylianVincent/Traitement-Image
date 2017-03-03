#include "pgm.h"
#include <math.h>

/*
    Detection de contours par filtre passe bande
    Utilisation de la TF du filtre LoG
 */

void filtrePasseBande(double **im1, double **im2, double **im3, double **im4, int nl, int nc, double sigma) {
    double factor;
    for(int u=0; u < nl; u++) {
        for(int v=0; v < nc; v++) {
            factor = (-4.0) * pow(M_PI, 2) * (pow(((u-(nl/2.0))/nl), 2) + pow(((v-(nc/2.0))/nc), 2)) * exp(-2 * pow(M_PI, 2) * pow(sigma, 2.0) * (pow(((u-(nl/2.0))/nl), 2) + pow(((v-(nc/2.0))/nc), 2)));
            im3[u][v] = im1[u][v] * factor;
            im4[u][v] = im2[u][v] * factor;
        }
    }
}

void passageParZero(double **im1, double **im2, int nl, int nc) {
    for(int u=0; u < nl-1; u++) {
        for(int v=0; v < nc-1; v++) {
            if (im1[u][v] * im1[u+1][v] < 0
                || im1[u][v] * im1[u][v+1] < 0) {
                im2[u][v] = 254;
            }
        }
    }
}

int main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
    int nb,nl,nc, oldnl,oldnc;
    unsigned char **im2=NULL,** im1=NULL;
    double** im4,** im5, ** im6, ** im7, **im8, **im9,**im10, **im11, **im12, **im13, **im14;

    if (ac < 3) {printf("Usage : %s entree sortie\n",av[0]); exit(1); }
    im1=lectureimagepgm(av[1],&nl,&nc);

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
    filtrePasseBande(im8, im9, im6, im7, nl, nc, 0.5);
    //filtrePasseBandeIdeal(im8, im9, im6, im7, nl, nc);
    fftshift(im6, im7, im8, im9, nl, nc);

    /* Calcul de la fft inverse de im8,im9 */
    ifft(im8,im9,im6,im7,nl,nc);

    passageParZero(im6, im4, nl, nc);

    // Sortie sur l'image im1
    /* Conversion en entier8bits de la partie reelle de la fftinverse,
    Suppresion des 0 qui ont servi a completer en utilisant la fonction crop
    Sauvegarde au format pgm de cette image qui doit etre identique a 'linverse video
    car on a realise la suite fftinv(fft(image))*/
    ecritureimagepgm(av[2],crop(imdouble2uchar(im4,nl,nc),0,0,oldnl,oldnc),oldnl,oldnc);
}
