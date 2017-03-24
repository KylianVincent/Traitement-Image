
#include "pgm.h"
#include "filters.h"

int main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
  int nb,nl,nc, oldnl,oldnc, sigma;
  unsigned char **im2=NULL,** im1=NULL;
  double** im4,** im5, ** im6, ** im7, **im8, **im9,**im10;

  if (ac < 3) {printf("Usage : %s entree sortie sigma\n",av[0]); exit(1); }
	/* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
  im1=lectureimagepgm(av[1],&nl,&nc);
  if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }

  sigma = av[2];
	/* transfomation en double -> calcul */ 
	im4 = imuchar2double(im1,nl,nc);
	/* Calcul des contours */
	bilateralFilter(im4,im5,nl,nc,1);

	/* transfomation en char -> sauvgade*/ 
	im2=imdouble2uchar(im4,nl,nc);
	/* Sauvegarde dans un fichier dont le nom est passe sur la ligne de commande */
  ecritureimagepgm(av[2],im2,nl,nc);
}
