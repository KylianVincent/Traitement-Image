#include "filters.h"
#include "math.h"
#include "float.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pgm.h"

/*
==========================================================
================ Filtrage non linéaire ===================
==========================================================
*/

/* ------------ Filtre médian ------------*/

int medianeFromHist(int n,int hist[]){
	int i,cummul = 0;
	for(i=0; i<256 && cummul<(n/2);i++){
		cummul+=hist[i];
	}
	return i;
}

void MedianFilter(double** sortie, double** entree, int nl, int nc, int n){
	//n : demi-taille du filtre
	//hist et la fréquence d'appaition de chaque valeur sur la zone filtré.
	//au début on calcule tout lhistogramme pour la première fenètre.
	//a chaque itération i,j on décale la fenetre, on retire les anciennes valeur et on ajoute les nouvelles.
	//pour trouver la median m c'est le plus petit m tel que : sum(i=0..m of hist[i])>n
	int l = 2*n+1;
	int hist[256]={0};
	//chaque ligne
	for(int i=0; i<nl; i++){
		//Initialisation de l'histogramme
		for(int k=0;k<256;k++)
			hist[k]=0;
		for(int k=0; k<l*l; k++){
			int in = (int)entree[prolongateByMirror(i-n+k/l,nl)][prolongateByMirror(k%l,nc)];
			hist[in]++;
			/*printf("%d\t%d\t%d\n",i-n+k/l,k%l,in);//*/
		}
		/*/	//	if(i==25)
			for(int k=0;k<256;k++)
				printf("%d\t%d\n",k,hist[k]);//*/
		sortie[i][0] = medianeFromHist(l*l, hist);
		//chaque colone
		for(int j=1; j<nc; j++){
			//MAJ de l'histogramme
			for(int k=0; k<l; k++){
				hist[(int)entree[prolongateByMirror(i-n+k,nl)][prolongateByMirror(j-n-1,nc)]]--;
				hist[(int)entree[prolongateByMirror(i-n+k,nl)][prolongateByMirror(j+n,nc)]]++;
				//printf("%d\t%d\t\t%d\t%d\n",i-n+k,j-n-1,i-n+k,j+n);//*/
			}
			/*printf("\n");//*/
			sortie[i][j] = medianeFromHist(l*l, hist);
			/*	if(j==52 && i==52){
			for(int k=0;k<256;k++)
				printf("%d\t%d\n",k,hist[k]);
			printf("%f\t",sortie[i][j]);
			}//*/
		}
	}
}

/* ------ Filtre adaptatif récursif ------*/
void adaptativeFilterInit(double** imSrc, double** imRes, double k, int nl, int nc) {
    int t = 0;
    //double oldDiff = DBL_MAX;

    //Terminaison condition
//    while (fabs(differenceBetweenImages(imSrc, imRes, nl, nc) - oldDiff) >= pow(10, -4)) {
    while ((fabs(differenceBetweenImages(imSrc, imRes, nl, nc)) >= 2) && t<201) {
        //printf("Diff : %f\n", fabs(differenceBetweenImages(imSrc, imRes, nl, nc)));
        t++;
        if (t%2 == 1) {
            adaptativeFilterRecursion(imSrc, imRes, k, nl, nc);
        } else {
            adaptativeFilterRecursion(imRes, imSrc, k, nl, nc);
        }
    }
        //printf("Diff : %f\n", fabs(differenceBetweenImages(imSrc, imRes, nl, nc)));

}

void adaptativeFilterRecursion(double** imSrc, double** imRes, double k, int nl, int nc) {
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

double median(int n, double x[]) {
	double temp;
	int i, j;
	// the following two loops sort the array x in ascending order
	for(i=0; i<n-1; i++) {
		for(j=i+1; j<n; j++) {
			if(x[j] < x[i]) {
				// swap elements
				temp = x[i];
				x[i] = x[j];
				x[j] = temp;
			}
		}
	}

	if(n%2==0) {
		// if there is an even number of elements, return mean of the two elements in the middle
		return((x[n/2] + x[n/2 - 1]) / 2.0);
	} else {
		// else return the element in the middle
		return x[n/2];
	}
}

void bilateralFilter(double** sortie, double** entree, int nl, int nc, double sigma1, double sigma2, bool med) {
	int x,y,i,j,ii,jj,k,l=(1+6*sigma1);
	double gaussienne[l*l+1],pixelNorme;
    if (sortie==NULL) sortie=alloue_image_double(nl,nc);
	//Pour chaque pixel du masque(sigma) -> précalcule de la gaussienne
	for(i=-3*sigma1; i<=3*sigma1; i++) {
		for(j=-3*sigma1; j<=3*sigma1; j++) {
			ii = (i+3*sigma1);
			jj = (j+3*sigma1);
			k = ii+l*jj;
			gaussienne[k] = exp((double)(i*i+j*j)/(-2*sigma1*sigma1));
			//printf("%f\t",gaussienne[k]);
		}
		//printf("\n");
	}
	//Pour tout les pixels de l'image
	for(x=0; x<nl; x++) {
		for(y=0; y<nc; y++) {
			pixelNorme=0;
			sortie[x][y] = 0;
			//Pour chaque pixel du masque(sigma)
			for(i=-3*sigma1; i<=3*sigma1; i++) {
				for(j=-3*sigma1; j<=3*sigma1; j++) {
					//printf("i%d,j%d,x%d,y%d\nx+i%d\n",i,j,x,y,x+i);
					ii = (i+3*sigma1);
					jj = (j+3*sigma1);
					k = ii+l*jj;
					double entreIJ = entree[prolongateByMirror(x+i,nl)][prolongateByMirror((y+j),nc)];
					double entreXY = entree[prolongateByMirror(x,nl)][prolongateByMirror((y),nc)];
					if(med) {
						double medIJ[] = {	entree[prolongateByMirror(x+i-1,nl)][prolongateByMirror(y+j-1,nc)],
											entree[prolongateByMirror(x+i-1,nl)][prolongateByMirror(y+j,nc)],
											entree[prolongateByMirror(x+i-1,nl)][prolongateByMirror(y+j+1,nc)],
											entree[prolongateByMirror(x+i,nl)][prolongateByMirror(y+j-1,nc)],
											entree[prolongateByMirror(x+i,nl)][prolongateByMirror(y+j,nc)],
											entree[prolongateByMirror(x+i,nl)][prolongateByMirror(y+j+1,nc)],
											entree[prolongateByMirror(x+i+1,nl)][prolongateByMirror(y+j-1,nc)],
											entree[prolongateByMirror(x+i+1,nl)][prolongateByMirror(y+j,nc)],
											entree[prolongateByMirror(x+i+1,nl)][prolongateByMirror(y+j+1,nc)], };
						entreIJ = median(9,medIJ);
						double medXY[] = {	entree[prolongateByMirror(x-1,nl)][prolongateByMirror(y-1,nc)],
											entree[prolongateByMirror(x-1,nl)][prolongateByMirror(y,nc)],
											entree[prolongateByMirror(x-1,nl)][prolongateByMirror(y+1,nc)],
											entree[prolongateByMirror(x,nl)][prolongateByMirror(y-1,nc)],
											entree[prolongateByMirror(x,nl)][prolongateByMirror(y,nc)],
											entree[prolongateByMirror(x,nl)][prolongateByMirror(y+1,nc)],
											entree[prolongateByMirror(x+1,nl)][prolongateByMirror(y-1,nc)],
											entree[prolongateByMirror(x+1,nl)][prolongateByMirror(y,nc)],
											entree[prolongateByMirror(x+1,nl)][prolongateByMirror(y+1,nc)], };
						entreXY = median(9,medXY);

					}
					double diff =entreIJ-entreXY;
					double intesitDiff = (diff*diff)/(double)(-2*sigma2*sigma2);
					double filtre = gaussienne[k]*exp(intesitDiff);
					pixelNorme += filtre;
					sortie[x][y] += entree[prolongateByMirror(x+i,nl)][prolongateByMirror((y+j),nc)]*filtre;
				}
			}
			//Normalisation
			sortie[x][y] /= pixelNorme;
		}
	}
}


/* ------------- Filtre patch ------------*/

void NLMeansFilter(double** imSrc, double** imRes, int nl, int nc, int t, int r, double sigma) {
    double w;
    double wPQ;
    double dPQ;
//    printf("T : %i\nR: %i\n", t, r);
    for (int u = 0; u < nl; u++) {
        for (int v = 0; v < nc; v++) {
            // Patch P centré en (u,v)
            imRes[u][v] = 0;
            wPQ = 0;
            for (int x = u-t; x < u+t; x++) {
                for(int y = v-t; y < v+t; y++) {
                    // Patch Q centré en (x,y)
                    dPQ = 0;
                    for (int i = -r; i <= r; i++) {
                        for (int j = -r; j <= r; j++) {
                            //printf(" Indices : (%i)%i - %i   ;    %i - %i\n", x+i, prolongateByMirror(x+i, nl), prolongateByMirror(y+j, nc), prolongateByMirror(u+i, nl), prolongateByMirror(v+j, nc));
                            dPQ += pow(imSrc[prolongateByMirror(x+i, nl)][prolongateByMirror(y+j, nc)] - imSrc[prolongateByMirror(u+i, nl)][prolongateByMirror(v+j, nc)], 2) * 1/pow(2*r+1, 2);
                        }
                    }
                    w = exp(-(dPQ/(2*pow(sigma, 2.0))));
                    wPQ += w;
                    //CALCULIMSORTIE
                    imRes[u][v] += w * imSrc[prolongateByMirror(x, nl)][prolongateByMirror(y, nc)];
                }
            }
            //CALCULIMSORTIE
            imRes[u][v] /= wPQ;
//            printf("%f\n", imRes[u][v]);
        }
    }
}

/* ---------- Extimation du bruit --------*/
double noiseEstimation(double ** im, int nl, int nc, int t, double p) {
    double** tmp=alloue_image_double(nl,nc);
    //Convolution par masque
    for (int u = 0; u < nl; u++) {
        for (int v = 0; v < nc; v++) {
            tmp[u][v] = im[u][v] - im[prolongateByMirror(u-1, nl)][v]
                        - im[prolongateByMirror(u+1, nl)][v]
                        - im[u][prolongateByMirror(v-1, nc)]
                        - im[u][prolongateByMirror(v+1, nc)];
        }
    }

    //Histogramme des variances
    int histVar[256*256]={0};
    int var;
    int moy;
    for (int u = 0; u < nl; u++) {
        for (int v = 0; v < nc; v++) {
            var = 0;
            moy = 0;
            for (int x = u-t; x <= u+t; x++) {
                for (int y = v-t; y <= v+t; y++) {
                    var += pow(tmp[prolongateByMirror(x, nl)][prolongateByMirror(y, nc)], 2);
                    moy += tmp[prolongateByMirror(x, nl)][prolongateByMirror(y, nc)];
                }
            }
            moy /= pow(2.0*t+1.0,2);
            var = var*(1.0/pow(2.0*t+1.0,2)) - pow(moy, 2);
            //printf("Var : %i\n", var);
            histVar[var]++;
        }
    }

    // for (int i = 0; i < 256*256; i++) {
    //     perror("%i : %i\n", i, histVar[i]);
    // }

    //Parcours de l'Histogramme
    int parcourus = 0;
    int aParcourir = pow(256, 2)*p/100;
    // Initialisation à -1 pour attaquer le premier tour de boucle à 0
    int varianceCour = -1;
    while (parcourus < aParcourir) {
        // printf("%i TO %i\n", parcourus, aParcourir);
        varianceCour++;
        parcourus += histVar[varianceCour];
    }
    // printf("%i TO %i\n", parcourus, aParcourir);
    return sqrt(1.13*varianceCour);
}

/* ---------------- Utils ----------------*/
int prolongateByMirror(int u, int nl) {
    if (u >= nl) {
        return nl - (u - nl + 1);
    } else if (u < 0) {
        return -(u + 1);
    }
    return u;
}

double differenceBetweenImages(double** im1, double** im2, int nl, int nc) {
    double diff = 0;
    for (int u = 0; u < nl; u++) {
        for (int v = 0; v < nc; v++) {
            diff += pow(im1[u][v] - im2[u][v], 2.0);
            //printf("%f\n", im1[u][v] - im2[u][v]);
        }
    }
    return sqrt(diff);
}
