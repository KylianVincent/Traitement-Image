/*
==========================================================
================ Filtrage non linéaire ===================
==========================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "pgm.h"

/* ------------ Filtre médian ------------*/

/* ------ Filtre adaptatif récursif ------*/

/* ----------- Filtre bilatéral ----------*/

void bilateralFilter(double** sortie, double** entree, int nl, int nc, int sigma1) { 
	int x,y,i,j,ii,jj,k,l=(1+6*sigma1),pixelNorme;
	double gaussienne[l*l+1];
    if (sortie==NULL) sortie=alloue_image_double(nl,nc);
	//Pour chaque pixel du masque(sigma) -> précalcule de la gaussienne
	for(i=-3*sigma1; i<3*sigma1; i++) {
		for(j=-3*sigma1; j<3*sigma1; j++) {
			ii = (i+3*sigma1);
			jj = (j+3*sigma1);
			k = ii+l*jj;
			gaussienne[k] = exp((i*i+j*j)/(-2*sigma1*sigma1));
		}
	}
	//Pour tout les pixels de l'image
    //for(y=0; y<nl; y++) {
      //  for(x=0; x<nc; x++) {
    for(y=30; y<nl-30; y++) {
        for(x=30; x<nc-30; x++) {
			pixelNorme=0;
			//Pour chaque pixel du masque(sigma)
			for(i=-3*sigma1; i<3*sigma1; i++) {
				for(j=-3*sigma1; j<3*sigma1; j++) {
					printf("i%d,j%d,x%d,y%d\nx+i%d\n",i,j,x,y,x+i);
					ii = (i+3*sigma1);
					jj = (j+3*sigma1);
					k = ii+l*jj;
					sortie[x][y] = gaussienne[k]*exp((i*i+j*j)/(-2*sigma1*sigma1));//TODO
					pixelNorme += sortie[x][y];
					sortie[x][y] *= entree[x+i][(y+j)];
				}
			}
			//Normalisation
			sortie[x][y] /= pixelNorme;
		}
	}
}


/* ------------- Filtre patch ------------*/

/* ---------- Extimation du bruit --------*/
