#include <stdbool.h>
/*
==========================================================
================ Filtrage non linéaire ===================
==========================================================
*/

/* ------ Filtre médian ------*/
void MedianFilter(double** sortie, double** entree, int nl, int nc, int n);

/* ------ Filtre adaptatif récursif ------*/
void adaptativeFilterInit(double** imSrc, double** imRes, double k, int nl, int nc);

void adaptativeFilterRecursion(double** imSrc, double** imRes, double k, int nl, int nc);

/* ------ Filtre bilatéral ------*/
void bilateralFilter(double** sortie, double** entree, int nl, int nc, double sigma1, double sigma2, bool med);

/* ------ Filtre patch ------*/
// t : Taille des régions
// r : Taille des patchs
void NLMeansFilter(double** imSrc, double** imRes, int nl, int nc, int t, int r, double sigma);

/* ------ Extimation du bruit ------*/
double noiseEstimation(double ** im, int nl, int nc, int t, double p);

/* ---------------- Utils ----------------*/
int prolongateByMirror(int u, int nl);

double differenceBetweenImages(double** im1, double** im2, int nl, int nc);
