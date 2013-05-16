/* lu.c
    Solution of the matrix equation  b = A x
    by LU (lower-upper -or- Cholesky) decomposition 
    and backsubstitution.
*/

#include <math.h>
#define TINY    (1.0e-20)

void ludcmp (double **a, int n, int *indx, double *d);
void lubksb (double **a, int n, int *indx, double *b);

/**************
 *** ludcmp ***
 **************
 * Cholesky decomposition of a matrix
 */
void ludcmp (double **a, int n, int *indx, double *d)
{
    int i, imax, j, k;
    double big, dum, sum, temp;
    double *vv;         /* stores the implicit doubling of each row */
    double  *vector (int nl, int nh);
    void    free_vector (double *v, int nl, int nh);
    void    nrerror (char error_text[]);

    vv = vector (1, n);
    *d = 1.0;           /* no row interchanges yet */
    /*
     * Loop over the row to get
     * implicit scaling information.
     */
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if ((temp=fabs(a[i][j])) > big)
                big = temp;
        if (big == 0.0)
            nrerror ("ludcmp: Singular matrix");
        /*
         * Otherwise, no non-zero largest element
         */
        vv[i] = 1.0/big;
    }
    /*
     * This is the loop over columns of Crout's method
     */
    for (j=1; j <= n; j++) {
        for (i = 1; i < j; i++) {   /* equation (2.3.12) except for i=j */
            sum = a[i][j];
            for (k = 1; k < i; k++)
                sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
        }
        /*
         * Initialize the search for the largest pivot element 
         */
        big = 0.0;
            /*
             * This is i=j of equation (2.3.12) and
             * i=j+1, ... N of equation (2.3.13)
             */
        for (i = j; i <= n; i++) {
            sum = a[i][j];
            for (k = 1; k < j; k++)
                sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            /*
             * Is the figure of merit for the pivot 
             * better than the best so far?
             */
             if ( (dum = vv[i]*fabs(sum)) >= big ) {
                big = dum;
                imax = i;
             }
        }
        if ( j != imax ) {      /* need to interchange rows? */
            for (k = 1; k <= n; k++) {  /* yes, do so ... */
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);         /* ... and interchange the parity of d */
            vv[imax] = vv[j];   /* ... and the scale factor */
        }
        indx[j] = imax;
        if (a[j][j] == 0.0)
            a[j][j] = TINY;
        /*
         * If the pivot element is sero, the matrix singular
         * (at least to the precision of the algorithm).
         * For some applications on singular matrices, it
         * is desirable to substitute TINY for zero.
         */
         if (j != n) {      /* finally, divide by the pivot element */
            dum = 1.0/(a[j][j]);
            for (i = j+1; i <= n; i++)
                a[i][j] *= dum;
         }
    }
    free_vector (vv, 1, n);
}


/**************
 *** lubksb ***
 **************
 * Solve the set of linear equations Ax = b.
 * "A" must be the LU decomposition of the matrix A.
 */
void lubksb (double **a, int n, int *indx, double *b)
{
    int i, ii=0, ip, j;
    double sum;

    /*
     * When ii is set to a positive value, it will become
     * the index of the first nonvanishing element of b.
     * We now do the forward substitution, equation (2.3.6).
     * The only new wrinkle is to unscramble the
     * permutation as we go.
     */
    for (i = 1; i <= n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii)
            for (j = ii; j <= i-1; j++)
                sum -= a[i][j]*b[j];
        else
            /*
             * A non-zero element was encountered,
             * so from now on we will have to do
             * the sums in the loop above.
             */
            if (sum)
                ii = i;
        b[i] = sum;
    }
    /*
     * Now we do the backsubstitution, equation (2.3.7)
     */
    for (i = n; i >= 1; i--) {
        sum = b[i];
        for (j = i+1; j <= n; j++)
            sum -= a[i][j]*b[j];
        /*
         * Store a component of the solution vector x 
         */
        b[i] = sum/a[i][i];
    }
}
