#include "pgm.h"
#include <math.h>
#include <stdbool.h>

        /*
                Calcul de l'inverse video d'une image. Si le parametre sortie est NULL, on cree une nouvelle image.
                On parcourt toute l'image ligen par ligne, colonne par colonne, on calcule son inverse video
                et on retourne l'image ainsi modifiee
        */
double** detection(double** sortie, double** entree, int nl, int nc) { 
	int i,j;
    if (sortie==NULL) sortie=alloue_image_double(nl,nc);
    for(i=1; i<nl-1; i++) {
        for(j=1; j<nc-1; j++) {
            sortie[i][j]=pow(entree[i+1][j+1]+entree[i][j+1]+entree[i-1][j+1]-(entree[i+1][j-1]+entree[i][j-1]+entree[i-1][j-1]),2);
            sortie[i][j]+=pow(entree[i+1][j-1]+entree[i+1][j]+entree[i+1][j+1]-(entree[i-1][j-1]+entree[i-1][j]+entree[i-1][j+1]),2);
            sortie[i][j]=sqrt(sortie[i][j]);
        }
    }
    return sortie;
}

double mean(double** entree, int nl, int nc) { 
	int i,j;
	double mean=0.0;
    for(i=0; i<nl; i++) {
        for(j=0; j<nc; j++) {
            mean+=entree[i][j];
		}
    }
	mean/=(i*j);
	return mean;
}

double** threshold(double threshold,double** sortie, double** entree, int nl, int nc) { 
	int i,j;
    if (sortie==NULL) sortie=alloue_image_double(nl,nc);
    for(i=0; i<nl; i++) {
        for(j=0; j<nc; j++) {
			if(entree[i][j] < threshold ) {
				sortie[i][j]=250.0;
			} else {
				sortie[i][j]=0.0;
			}
        }
    }
    return sortie;
}

bool isLocalMin(double** entree, int i, int j) {
	double threshold =1.0;
	return	entree[i][j] > entree[i-1][j]+threshold && entree[i][j] > entree[i+1][j]+threshold 
		|| 	entree[i][j] > entree[i][j-1]+threshold && entree[i][j] > entree[i][j+1]+threshold 
		|| 	entree[i][j] > entree[i-1][j+1]+threshold && entree[i][j] > entree[i+1][j-1]+threshold 
		|| 	entree[i][j] > entree[i-1][j-1]+threshold && entree[i][j] > entree[i+1][j+1]+threshold ;
}

double** onlyLocalMin(double** sortie, double** entree, int nl, int nc) { 
	int i,j;
    if (sortie==NULL) sortie=alloue_image_double(nl,nc);
    for(i=1; i<nl-1; i++) {
        for(j=1; j<nc-1; j++) {
			if(!isLocalMin(entree,i,j)) {
				sortie[i][j]=250.0;
			} else {
				sortie[i][j]=0.0;
			}
        }
    }
    return sortie;
}
	/*
		Exemple de code avec Entrees Sortie et transformations simples d'images
		S'utilise sous la forme  "exemple tangram.pgm res.pgm"
 	*/
int main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
  int nb,nl,nc, oldnl,oldnc;
  unsigned char **im2=NULL,** im1=NULL;
  double** im4,** im5, ** im6, ** im7, **im8, **im9,**im10;

  if (ac < 2) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
	/* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
  im1=lectureimagepgm(av[1],&nl,&nc);
  if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }

	/* transfomation en double -> calcul */ 
	im4 = imuchar2double(im1,nl,nc);
	/* Calcul des contours */
  im5=detection(NULL,im4,nl,nc);
  im4=threshold(mean(im5,nl,nc),NULL,im5,nl,nc);
  //im4=threshold(20.0,NULL,im5,nl,nc);
  //im4=onlyLocalMin(NULL,im5,nl,nc);

	/* transfomation en char -> sauvgade*/ 
	im2=imdouble2uchar(im4,nl,nc);
	/* Sauvegarde dans un fichier dont le nom est passe sur la ligne de commande */
  ecritureimagepgm(av[2],im2,nl,nc);
}
