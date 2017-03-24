/*
==========================================================
====================== Débruitage ========================
==========================================================
*/

void filtrageFrequentiel(double** imReal, double** imImag, double sigma, int nl, int nc);

void filtreSeparable(double** imSrc, double** imRes, double sigma, int nl, int nc, int w, double* tabCoef, int dirX);

void filtreBrut(double** imSrc, double** imRes, double sigma, int nl, int nc, int w, double* tabCoef);

void filtrageSpatial(double** imSrc, double** imRes, double sigma, int nl, int nc, int w);



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

void filtrePasseBandeFrequentiel(double **imReal, double **imImag, int nl, int nc, double sigma);

void passageParZero(double **imSrc, double **imRes, int nl, int nc);
