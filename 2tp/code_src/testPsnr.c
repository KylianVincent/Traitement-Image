#include "pgm.h"
#include "divers.c"

  /*
		Calcul du psnr sur deux images données en entrée
		S'utilise sous la forme  "exemple source.pgm generee.pgm"
 	*/
int main(int ac, char **av){

    int nl,nc;

    if (ac < 3) {printf("Usage : %s source generee\n",av[0]); exit(1); }
	   /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
    unsigned char** im1=NULL,** im2=NULL;
    im1=lectureimagepgm(av[1],&nl,&nc);
    im2=lectureimagepgm(av[2],&nl,&nc);

    double**im3=imuchar2double(im1,nl,nc);
    double**im4=imuchar2double(im2,nl,nc);

    if (im1==NULL || im2==NULL)  { puts("Lecture image impossible"); exit(1); }

    double psnrValue;
    psnrValue =  psnr_double(im3, im4, nl, nc);

    printf("%f\n", psnrValue);
}
