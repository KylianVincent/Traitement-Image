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
	double gaussienne[l*l];
    if (sortie==NULL) sortie=alloue_image_double(nl,nc);
	//Pour chaque pixel du masque(sigma) -> précalcule de la gaussienne
	for(i=-3*sigma1; i<3*sigma1; i++) {
		for(j=-3*sigma1; j<3*sigma1; j++) {
			ii = (i+3*sigma1);
			jj = (j+3*sigma1);
			k = ii+l*jj;
			//gaussienne[k] = exp((i*i+j*j)/(-2*sigma1*sigma1));
		}
	}
	//Pour tout les pixels de l'image
    for(y=6*sigma1+1; y<nl-6*sigma1-1; y++) {
        for(x=6*sigma1+1; x<nc-6*sigma1-1; x++) {
			pixelNorme=0;
			//Pour chaque pixel du masque(sigma)
			for(i=-3*sigma1; i<3*sigma1; i++) {
				for(j=-3*sigma1; j<3*sigma1; j++) {
					ii = (i+3*sigma1);
					jj = (j+3*sigma1);
					k = ii+l*jj;
					sortie[x][y] = gaussienne[k]*exp((i*i+j*j)/(-2*sigma1*sigma1));
					pixelNorme += sortie[x][y];
					sortie[x][y] *= entree[x+i][y+j];
				}
			}
			//Normalisation
			sortie[x][y] /= pixelNorme;
		}
	}
}


/* ------------- Filtre patch ------------*/

/* ---------- Extimation du bruit --------*/
