/* recipes.h

    Routines from Numerical Recipes
*/

#ifndef _RECIPES_
#define _RECIPES_

double  *vector(int nl, int nh);
void    free_vector (double *v, int nl, int nh);
int     *ivector(int nl, int nh);
void    free_ivector (int *v, int nl, int nh);
double  **matrix(int nrl, int nrh, int ncl, int nch);
void    free_matrix (double **m, int nrl, int nrh, int ncl, int nch);
void    hunt(double xx[], int n, double x, int *jlo);
void    nrerror (char error_text[]);
void    svdcmp(double **a, int m, int n, double w[], double **v);
void    svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
double  pythag(double a, double b);
double  gasdev(long *idum);
double  ran1(long *idum);

#define MAX(_A_,_B_) ((_A_)>(_B_)?(_A_):(_B_))
#define MIN(_A_,_B_) ((_A_)<(_B_)?(_A_):(_B_))
#define SQR(a) (pow((a),2))  /* let the compiler optimize this */
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#endif
