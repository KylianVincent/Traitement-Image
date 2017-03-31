/*
==========================================================
================ Filtrage non linéaire ===================
==========================================================
*/

/* ------ Filtre médian ------*/

/* ------ Filtre adaptatif récursif ------*/
void adaptativeFilterInit(double** imSrc, double** imRes, double k, int nl, int nc);

void adaptativeFilterRecursion(double** imSrc, double** imRes, double k, int nl, int nc);

/* ------ Filtre bilatéral ------*/
void bilateralFilter(double** sortie, double** entree, int nl, int nc, int sigma1);

/* ------ Filtre patch ------*/

/* ------ Extimation du bruit ------*/

/* ---------------- Utils ----------------*/
int prolongateByMirror(int u, int nl);

double differenceBetweenImages(double** im1, double** im2, int nl, int nc);
