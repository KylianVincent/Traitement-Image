#include "pgm.h"
#include <math.h>

void filtreGaussien(double** im1, double** im2, double** im3, double** im4, double sigma, int nl, int nc) {
    if (*im3==NULL) im3=imuchar2double(alloue_image(nl,nc), nl, nc);
    if (*im4==NULL) im4=imuchar2double(alloue_image(nl,nc), nl, nc);
    double gauss = -2.0*(pow(M_PI, 2.0)*(pow(sigma, 2)));
    for(int u=0; u<nl; u++) {
        for(int v=0; v<nc; v++) {
            im3[u][v] = im1[u][v] * exp(gauss * (pow(((u-(nl/2.0))/nl), 2) + pow(((v-(nc/2.0))/nc), 2)));
            im4[u][v] = im2[u][v] * exp(gauss * (pow(((u-(nl/2.0))/nl), 2) + pow(((v-(nc/2.0))/nc), 2)));
        }
    }
}

	/*
		Debruitage par la methode de gauss
		S'utilise sous la forme  "exemple tangram.pgm res.pgm sigma"
 	*/
int main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
  int nb,nl,nc, oldnl,oldnc;
  unsigned char **im2=NULL,** im1=NULL;
  double** im4,** im5, ** im6, ** im7, **im8, **im9,**im10, **im11, **im12, **im13, **im14;

  if (ac < 4) {printf("Usage : %s entree sortie sigma \n",av[0]); exit(1); }
	/* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
  im1=lectureimagepgm(av[1],&nl,&nc);
  const double SIGMA = atof(av[3]);
  if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }
	/* Calcul de son inverse video */
  double**im3=imuchar2double(im1,nl,nc);
  oldnl=nl; oldnc=nc;
	/*  la fft demande des puissances de 2. On padde avec des 0, mais les dimensions nl et nc changent */
  im7=padimdforfft(im3,&nl,&nc);
  /*
	On peut aussi directement utiliser
	im7=padimucforfft(im1,&nl,&nc);
	sans convertir im1 en image de réels
  */
	/* Creation des images pour les parties reelles et imagianires des fft  */
  im4=alloue_image_double(nl,nc); im5=alloue_image_double(nl,nc); im6=alloue_image_double(nl,nc);
	/* Calcul de la fft de im7,im4 */
  fft(im7,im4,im5,im6,nl,nc);
  /* Creation des images pour les parties reelles et imagianires des fft inverses */

  im9=alloue_image_double(nl,nc); im10=alloue_image_double(nl,nc); im11=alloue_image_double(nl,nc); im12=alloue_image_double(nl,nc); im13=alloue_image_double(nl,nc); im14=alloue_image_double(nl,nc);
  /* Application du filtre gaussien */
  fftshift(im5, im6, im9, im10, nl, nc);
  filtreGaussien(im9, im10, im11, im12, SIGMA, nl, nc);
  //printf("%f\n", im11[30][30]);
  fftshift(im11, im12, im13, im14, nl, nc);

        /* Calcul de la fft inverse de im13,im14 */
  ifft(im13,im14,im11,im12,nl,nc);
	/* Conversion en entier8bits de la partie reelle de la fftinverse,
	   Suppresion des 0 qui ont servi a completer en utilisant la fonction crop
	   Sauvegarde au format pgm de cette image qui doit etre identique a 'linverse video
	   car on a realise la suite fftinv(fft(image))*/
  ecritureimagepgm(av[2],crop(imdouble2uchar(im11,nl,nc),0,0,oldnl,oldnc),oldnl,oldnc);
}
