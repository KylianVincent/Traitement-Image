#include "filters.h"
#include "math.h"
#include <stdbool.h>
#include <stdio.h>

/*
==========================================================
================ Filtrage non linéaire ===================
==========================================================
*/

/* ------------ Filtre médian ------------*/

/* ------ Filtre adaptatif récursif ------*/
void adaptativeFilterInit(double** imSrc, double** imRes, double k, int nl, int nc) {
    int t = 0;

    //Terminaison condition
    while (t < 30) {
        t++;
        if (t%2 == 1) {
            adaptativeFilterRecursion(imSrc, imRes, k, nl, nc, t);
        } else {
            adaptativeFilterRecursion(imRes, imSrc, k, nl, nc, t);
        }
    }
}

void adaptativeFilterRecursion(double** imSrc, double** imRes, double k, int nl, int nc, int t) {
    int um1, up1, vm1, vp1;
    double omega[nl][nc];
    double num, denom;

    //Omega_t computation
    //For each pixel
    for (int u = 0; u<nl; u++) {
        for (int v = 0; v<nc; v++) {
            //Mirror-prolongated image
            um1 = prolongateByMirror(u-1, nl);
            up1 = prolongateByMirror(u+1, nl);
            vm1 = prolongateByMirror(v-1, nc);
            vp1 = prolongateByMirror(v+1, nc);

            //printf("%d - %d - %d - %d\n", um1, up1, vm1, vp1);

            omega[u][v]=exp(-(pow(imSrc[up1][v] - imSrc[um1][v], 2.0) + pow(imSrc[u][vp1] - imSrc[u][vm1], 2.0))/(2*pow(k, 2.0)));
        }
    }

    //Image computation
    //For each pixel
    for (int u = 0; u<nl; u++) {
        for (int v = 0; v<nc; v++) {
            num = 0;
            denom = 0;
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    num += omega[prolongateByMirror(u+i, nl)][prolongateByMirror(v+j, nc)] * imSrc[prolongateByMirror(u+i, nl)][prolongateByMirror(v+j, nc)];
                    denom += omega[prolongateByMirror(u+i, nl)][prolongateByMirror(v+j, nc)];
                }
            }
            imRes[u][v] = num/denom;
        }
    }
}

/* ----------- Filtre bilatéral ----------*/

//void bilateralFilter(double** sortie, double** entree, int nl, int nc, int sigma1) {
//	int x,y,i,j,pixelNorm;
//	double gaussienne[pow((1+6*sigma1),2)];
//    if (sortie==NULL) sortie=alloue_image_double(nl,nc);
//	//Pour chaque pixel du masque(sigma) -> précalcule de la gaussienne
//	for(i=-3*sigma1; i<3*sigma1; i++) {
//		for(j=-3*sigma1; j<3*sigma1; j++) {
//			//TODO calcul gaussienne
//		}
//	}
//	//Pour tout les pixels de l'image
//    for(y=1; y<nl-1; y++) {
//        for(x=1; x<nc-1; x++) {
//			pixelNorm=0;
//			//Pour chaque pixel du masque(sigma)
//			for(i=-3*sigma1; i<3*sigma1; i++) {
//				for(j=-3*sigma1; j<3*sigma1; j++) {
//					//TODO calcul sortit + pixelNorme
//					/*
//            sortie[i][j]=pow(entree[i+1][j+1]+entree[i][j+1]+entree[i-1][j+1]-(entree[i+1][j-1]+entree[i][j-1]+entree[i-1][j-1]),2);
//            sortie[i][j]+=pow(entree[i+1][j-1]+entree[i+1][j]+entree[i+1][j+1]-(entree[i-1][j-1]+entree[i-1][j]+entree[i-1][j+1]),2);
//            sortie[i][j]=sqrt(sortie[i][j]);
//			*/
//				}
//			}
//			//Normalisation
//			sortie[x][y] /= pixelNorm;
//		}
//	}
//}


/* ------------- Filtre patch ------------*/

/* ---------- Extimation du bruit --------*/

/* ---------------- Utils ----------------*/
int prolongateByMirror(int u, int nl) {
    if (u >= nl) {
        return nl - (nl - u + 1);
    } else if (u < 0) {
        return -(u + 1);
    }
    return u;
}
