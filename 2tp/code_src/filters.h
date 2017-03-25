/*
==========================================================
================ Filtrage non linéaire ===================
==========================================================
*/

/* ------ Filtre médian ------*/

/* ------ Filtre adaptatif récursif ------*/
void adaptativeFilterInit(double** imSrc, double** imRes, double k, int nl, int nc);

void adaptativeFilterRecursion(double** imSrc, double** imRes, double k, int nl, int nc, int t);

/* ------ Filtre bilatéral ------*/

/* ------ Filtre patch ------*/

/* ------ Extimation du bruit ------*/

/* ---------------- Utils ----------------*/
int prolongateByMirror(int u, int nl);
