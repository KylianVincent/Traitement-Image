#include <stdbool.h>
/*
==========================================================
================ Filtrage non linéaire ===================
==========================================================
*/

/* ------ Filtre médian ------*/
void MedianFilter(double** sortie, double** entree, int nl, int nc);

/* ------ Filtre adaptatif récursif ------*/
void adaptativeFilterInit(double** imSrc, double** imRes, double k, int nl, int nc);

void adaptativeFilterRecursion(double** imSrc, double** imRes, double k, int nl, int nc, int t);

/* ------ Filtre bilatéral ------*/
void bilateralFilter(double** sortie, double** entree, int nl, int nc, double sigma1, double sigma2, bool med);

/* ------ Filtre patch ------*/

/* ------ Extimation du bruit ------*/

/* ---------------- Utils ----------------*/
int prolongateByMirror(int u, int nl);
