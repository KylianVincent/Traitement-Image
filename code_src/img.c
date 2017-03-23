#include "img.h"
#include "pgm.h"
#include <math.h>
/*
==========================================================
====================== Débruitage ========================
==========================================================
*/

void filtrageFrequentiel(double** imReal, double** imImag, double sigma, int nl, int nc) {
    double gauss = -2.0*(pow(M_PI, 2.0)*(pow(sigma, 2)));
    for(int u=0; u<nl; u++) {
        for(int v=0; v<nc; v++) {
            imReal[u][v] = imReal[u][v] * exp(gauss * (pow(((u-(nl/2.0))/nl), 2) + pow(((v-(nc/2.0))/nc), 2)));
            imImag[u][v] = imImag[u][v] * exp(gauss * (pow(((u-(nl/2.0))/nl), 2) + pow(((v-(nc/2.0))/nc), 2)));
        }
    }
}

void filtreSeparable(double** imSrc, double** imRes, double sigma, int nl, int nc, int w, double* tabCoef, int dirX) {
    //Filtrage spatial
    for(int u=0; u < nl; u++) {
        for(int v=0; v < nc; v++) {
            imRes[u][v] = 0.0;
            for (int j = -w; j < w; j++) {
                // On prolonge l'image par périodicité
                imRes[u][v] += tabCoef[abs(j)] * imSrc[(u + (dirX == 1 ? j:0) + nl) % nl][(v + (dirX == 1 ? 0:j) + nc) % nc];
            }
        }
    }
}

void filtreBrut(double** imSrc, double** imRes, double sigma, int nl, int nc, int w, double* tabCoef) {
    //Filtrage spatial
    for(int u=0; u<nl; u++) {
        for(int v=0; v<nc; v++) {
            for (int j = -w; j<w; j++) {
                for (int i = -w; i<w; i++) {
                    // On prolonge l'image par périodicité
                    imRes[u][v] += tabCoef[abs(j)] * tabCoef[abs(i)] * imSrc[(u + i + nl) % nl][(v + j + nc) % nc];
                }
            }
        }
    }
}

void filtrageSpatial(double** imSrc, double** imRes, double sigma, int nl, int nc, int w) {
    //Precalcul des coefficients exponentiels
    double tabCoef[w];
    for (int i = 0; i<w; i++) {
        tabCoef[i] = exp(-(pow(i, 2.0) / (2.0 * pow(sigma, 2))));
    }
    double** tmp=alloue_image_double(nl,nc);

    //Filtrage utilisant la proprriété de séparation
    filtreSeparable(imSrc, tmp, sigma, nl, nc, w, tabCoef, 1);
    filtreSeparable(tmp, imRes, sigma, nl, nc, w, tabCoef, 0);
    //Filtrage ne l'utilisant pas
//    filtreBrut(imSrc, imRes, sigma, nl, nc, w, tabCoef);

    //Multiplication par le facteur
    for(int u=0; u < nl; u++) {
        for (int v = 0; v < nc; v++) {
            imRes[u][v] *= 1 / (2 * M_PI * pow(sigma, 2.0));
        }
    }
}



/*
==========================================================
======================== Contours ========================
==========================================================
*/
/*
============== Gradiant method =================
*/








/*
============== Laplacian Filter ================
*/

void filtrePasseBandeFrequentiel(double **imReal, double **imImag, int nl, int nc, double sigma) {
    double factor;
    for(int u=0; u < nl; u++) {
        for(int v=0; v < nc; v++) {
            factor = (-4.0) * pow(M_PI, 2) * (pow(((u-(nl/2.0))/nl), 2) + pow(((v-(nc/2.0))/nc), 2)) * exp(-2 * pow(M_PI, 2) * pow(sigma, 2.0) * (pow(((u-(nl/2.0))/nl), 2) + pow(((v-(nc/2.0))/nc), 2)));
            imReal[u][v] = imReal[u][v] * factor;
            imImag[u][v] = imImag[u][v] * factor;
        }
    }
}

void passageParZero(double **imSrc, double **imRes, int nl, int nc) {
    for(int u=0; u < nl-1; u++) {
        for(int v=0; v < nc-1; v++) {
            if ((imSrc[u][v] * imSrc[u+1][v] <= -0 && abs(imSrc[u][v] - imSrc[u+1][v]) >= 1)
                || (imSrc[u][v] * imSrc[u][v+1] <= -0 && abs(imSrc[u][v] - imSrc[u][v+1]) >= 1)) {
                imRes[u][v] = 0;
            } else {
                imRes[u][v] = 254;
            }
        }
    }
}
