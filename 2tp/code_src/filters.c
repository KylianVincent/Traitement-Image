/*
==========================================================
================ Filtrage non linéaire ===================
==========================================================
*/
#include <math.h>
#include <stdbool.h>

/* ------------ Filtre médian ------------*/

/* ------ Filtre adaptatif récursif ------*/

/* ----------- Filtre bilatéral ----------*/

void bilateralFilter(double** sortie, double** entree, int nl, int nc, int sigma1) { 
	int x,y,i,j,pixelNorm;
	double gaussienne[pow((1+6*sigma1),2)];
    if (sortie==NULL) sortie=alloue_image_double(nl,nc);
	//Pour chaque pixel du masque(sigma) -> précalcule de la gaussienne
	for(i=-3*sigma1; i<3*sigma1; i++) {
		for(j=-3*sigma1; j<3*sigma1; j++) {
			//TODO calcul gaussienne
		}
	}
	//Pour tout les pixels de l'image
    for(y=1; y<nl-1; y++) {
        for(x=1; x<nc-1; x++) {
			pixelNorm=0;
			//Pour chaque pixel du masque(sigma)
			for(i=-3*sigma1; i<3*sigma1; i++) {
				for(j=-3*sigma1; j<3*sigma1; j++) {
					//TODO calcul sortit + pixelNorme
					/*
            sortie[i][j]=pow(entree[i+1][j+1]+entree[i][j+1]+entree[i-1][j+1]-(entree[i+1][j-1]+entree[i][j-1]+entree[i-1][j-1]),2);
            sortie[i][j]+=pow(entree[i+1][j-1]+entree[i+1][j]+entree[i+1][j+1]-(entree[i-1][j-1]+entree[i-1][j]+entree[i-1][j+1]),2);
            sortie[i][j]=sqrt(sortie[i][j]);
			*/
				}
			}
			//Normalisation
			sortie[x][y] /= pixelNorm;
		}
	}
}


/* ------------- Filtre patch ------------*/

/* ---------- Extimation du bruit --------*/
