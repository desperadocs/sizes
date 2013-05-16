/* regular.c

Regularization solution of size distributions from SAS data

Solves the linear equation
	I[1..M] +/- s[1..M] = G[1..M][1..N] x[1..N]
by minimization of the functional
	Q = a S + C
where the fit to the data (chiSquared) is given by
	C = || I - G x ||
and the applied constraint is
	S = | x |		or
	S = | x' |		or
	S = | x'' |		or ...
subject to the constraint that
	C = M +/- e * sqrt(2M)
where e ~ 0.1 or so.

Minimization of Q leads to the equation
	B[1..N] = A[1..N][1..N] x[1..N]
which is linear in x[].  Here,
	B[i] = Sum[ G[i][j] * I[j] / (s[j]^2), {j,1,M}]
The matrix A[][] is constructed with different values of "a" as
	A[1..N][1..N] = D[1..N][1..N] + a * H[1..N][1..N]
where
	D[i][k] = Sum[ G[i][j] * G[k][j] / (s[j]^2), {j,1,M}]]
and H[i][k] is dependent on the applied constraint.
The preferred constraint
	S = | x'' | = Sum[ (2x[i] - x[i-1] - x[i+1])^2, {i, 2, N-2}]
minimizes the curvature of the solution (as does cubic splines).
Then,
H[1..N][1..N] = 
	|  1  -2   1   0   0   0   0   0            ...        0 |
	| -2   5  -4   1   0   0   0   0            ...        0 |
	|  1  -4   6  -4   1   0   0   0            ...        0 |
	|  0   1  -4   6  -4   1   0   0            ...        0 |
	|  0   0   1  -4   6  -4   1   0            ...        0 |
	|                  ...                                   |
	|  0               ...         0   1  -4   6  -4   1   0 |
	|  0               ...         0   0   1  -4   6  -4   1 |
	|  0               ...         0   0   0   1  -4   5  -2 |
	|  0               ...         0   0   0   0   1  -2   1 |
which is a symmetric, banded matrix.


The solution algorithm is as follows:
	 1. Measure I[] +/- s[]
	 2. assume a binning arrangement for x[]
	 3. Calculate G[][] from given description of the system
	 4. Calculate B[], D[][], and H[][] from G[][], I[], and s[]
	 5. assume a trial "a"
	 6. calculate A[][] = D[][] + a * H[][]
	 7. LU decompose A[][]
	 8. Backsubstitute B[] to get x[] = Inverse(A[][]) B[]
	 9. Calculate C
	10. C == M +/- e * sqrt(2M) ?
			no:  repeat from step 5 with better "a"
			yes: accept solution of x[]
	11. Calculate dx[1..N] from covariance matrix

Comments:
	Steps 5-10 can be incorporated into a search for "a"
	by bisection.  It is known that -80 < log (a) < +80.
	When "a" is too large, then C will be large, indicating
	that smoothing is turned ON.  When "a" is too low, C will
	be small, indicating the least squares solution.
	The bisection search goes as follows:
		upper = 80.0;
		lower = -80.0;
		tolerance = e * sqrt(2*M);
		do {
			midpoint = (upper + lower) * 0.5;
			a = pow (10.0, midpoint);
			steps 6-9
			if (C > M)
				upper = midpoint;
			else
				lower = midpoint;
		} while ( fabs (C-M) > tolerance );  // step 10

	Steps 7 and 8 can be replaced by singular value decomposition
	and backsubstitution, respectively.  However, this is unnecessary
	because the "a" term serves the same purpose as the zeroing
	of the singular values.  Thus it takes too long.
	
	Lawson and Hanson have a routine that will solve
		B[] = A[][] x[]
	while constraining x[i] >= 0.0 for i=[1..N].  This is routine
	NNLS in their book, "Solving Least Squares Problems,"
	(c)1971, Prentice-Hall, New York.

********************************************************************/

#include <stdio.h>
#include <math.h>
#include "nnls2.h"

#define EPSILON	(0.1)
/* #define SCAN_ALPHA_RANGE	/* for scanning across the interesting part */
#define UPPER_ALPHA		(80.0)
#define LOWER_ALPHA		(-UPPER_ALPHA)
#define STEP_ALPHA		(1)
#define ITERATION_MAXIMUM       (32)
#define REG_NNLS_METHOD          2
#define REG_LU_METHOD            0
#define MAXENT_METHOD            1

void	ludcmp (double **a, int n, int *indx, double *d);
void	lubksb (double **a, int n, int *indx, double *b);
double  *vector(int nl, int nh);
void    free_vector (double *v, int nl, int nh);
int     *ivector(int nl, int nh);
void    free_ivector (int *v, int nl, int nh);
double  **matrix(int nrl, int nrh, int ncl, int nch);
void    free_matrix (double **m, int nrl, int nrh, int ncl, int nch);
void    call_nnls(double **aa, int m, int n, double *b, double *x);

int Regularize_v2 (double *I, double *s, int M, 
	           double *r, double *x, int N, double **G, int maxIter,
                   int method
    );

/*******************************************
 ******************************************* Solve
 *******************************************
 *
 * Find the size distribution x(r), given the particular
 * value of the LaGrange multiplier, 10^midpoint
 */
static void Solve (
	double midpoint,
	double *I, double *s, int M, 
	double *r, double *x, int N, double **G,
	double *ChiSqr, double *S,
	double *B, double **A, double **D, double **H, double *w, double *d,
	int *indx, double *a,
        int method
)
{
	double temp;
	int i, j, k;
	*a = pow ((double) 10.0, (double) midpoint);
	for (i = 1; i <= N; i++) {
		for (k = 1; k <= N; k++)
			A[i][k] = D[i][k]  +  (*a) * H[i][k];		/* step 6 */
	}
        switch (method) {
          case REG_NNLS_METHOD:
	    /* Non-Negative Least Squares solution of Ax=B */
            call_nnls((double **) A, N, N, B, x);
            break;
          default:
          case REG_LU_METHOD:
	    /* LU decomposition and backsubstitution solution of Ax=B */
            ludcmp (A, N, indx, &temp);     /* decomposition of A[][] */
	    for (i = 1; i <= N; i++) x[i] = B[i];   /* copy for in-place solution */
	    lubksb (A, N, indx, x);         /* backsubstitution with B[] --> x[] */
            break;
        }
	*ChiSqr = 0.0;
	for (j = 1; j <= M; j++) {	/* calculate the misfit */
		temp = -I[j];
		for (i = 1; i <= N; i++)
			temp += G[i][j] * x[i];		/* the equation solved */
		*ChiSqr += temp*temp * w[j];
	}
	*S = 0.0;
	for (i = 2; i <= N-1; i++) {
		temp = 2*x[i] - x[i-1] - x[i+1];
		*S += temp*temp;
	}
}


/*******************************************
 ******************************************* Regularize_v2
 *******************************************
 *
 * Perform the regularization treatment, as described above
 */
int Regularize_v2 (
     double *I, double *s, int M,    /* measured data */
     double *r, double *x, int N,    /* solution vector */
     double **G,                     /* design matrix */
     int maxIter,                    /* give up after this many tries */
     int method)                     /* method=0: LU  method=2: NNLS */
{
	int i, j, k;
	double *B, **A, **D, **H, *w, *d;
	double a, b, upper, lower, midpoint;
	double tolerance, temp, ChiSqr, S;
	int *indx;
	int iterationMaximum;

		/*
		 * allocate the memory to solve the problem
		 */
	indx = ivector (1, N);
	B = vector (1, N);
	w = vector (1, M);
	A = matrix (1, N, 1, N);
	D = matrix (1, N, 1, N);
	H = matrix (1, N, 1, N);
		/*
		 * perform the initial calculations
		 */
	for (j = 1; j <= M; j++)
		w[j] = 1.0 / (s[j]*s[j]);	/* weights are inverse variances */
	for (i = 1; i <= N; i++) {
		d = &B[i];				/* for addressing speed, use a pointer */
		*d = 0.0;
		for (j = 1; j <= M; j++)
			*d += G[i][j] * I[j] * w[j];
		for (k = 1; k <= N; k++) {
			d = &(D[i][k]);
			*d = 0.0;
			for (j = 1; j <= M; j++)
				*d += G[i][j] * G[k][j] * w[j];
			H[i][k] = 0.0;
		}
		if ( i > 2 && i < N-1 ) {	/* this constrains curvature, x'' */
			H[i][i-2] =  1.0;
			H[i][i-1] = -4.0;
			H[i][i]   =  6.0;
			H[i][i+1] = -4.0;
			H[i][i+2] =  1.0;
		}
	}
		/* special cases */
	H[1][1] = 1.0;	H[1][2] = -2.0;    H[1][3] = 1.0;
	H[N][N] = 1.0;	H[N][N-1] = -2.0;  H[N][N-2] = 1.0;
	H[2][1] = -2.0; 
		H[2][2] = 5.0;  
			H[2][3] = -4.0;  
				H[2][4] = 1.0;
	H[N-1][N] = -2.0; 
		H[N-1][N-1] = 5.0;  
			H[N-1][N-2] = -4.0;  
				H[N-1][N-3] = 1.0;
		/* also, constrain the magnitude of the first & last points
		 * this, in effect, forces the distribution to zero at the ends
		 */
	H[1][1] += 1.0;
	H[N][N] += 1.0;
		/*
		 * search for the optimal "a" by bisection
		 */
	tolerance = EPSILON * sqrt (2.0 * (double) M);
	upper = UPPER_ALPHA;		/* upper bound for search */
	lower = LOWER_ALPHA;		/* lower bound for search */
	printf ("%15s    %15s    %15s    %15s\n", 
			"log10(a)", "ChiSqr", "S", "Q=ChiSqr + a*S");
        fflush (stdout);

#ifdef SCAN_ALPHA_RANGE
	for (midpoint = LOWER_ALPHA; 
			midpoint <= UPPER_ALPHA; 
			midpoint += STEP_ALPHA) {
		char msg[256];
		FILE *path;

		Solve (midpoint, I,s,M, r,x,N, G, &ChiSqr, 
                       &S, B,A,D,H,w,d,indx, &a, method);

#ifdef __MWERKS__
		printf ("%15g    %15g    %15g    %15g\n",
			   midpoint, ChiSqr,  S,   ChiSqr + a*S);
		sprintf (msg, "f(%g)", midpoint);
		path = fopen (msg, "w");
		if (path) {
			fprintf (path, "%g\t : log10(alpha)\n", midpoint);
			fprintf (path, "%g\t : ChiSqr\n", ChiSqr);
			fprintf (path, "%g\t : S\n", S);
			fprintf (path, "%g\t : ChiSqr + a*S\n", ChiSqr + a*S);
			fprintf (path, "%s\t%s\n", "r", "f*dr");
			for (i = 1; i <= N; i++)
				fprintf (path, "%g\t%g\n", r[i], x[i]);
			fclose (path);
		}
#else
		printf ("%15lg    %15lg    %15lg    %15lg\n",
			   midpoint,  ChiSqr,   S,    ChiSqr + a*S);
		sprintf (msg, "f(%lg)", midpoint);
		path = fopen (msg, "w");
		if (path) {
			fprintf (path, "%lg\t : log10(alpha)\n", midpoint);
			fprintf (path, "%lg\t : ChiSqr\n", ChiSqr);
			fprintf (path, "%lg\t : S\n", S);
			fprintf (path, "%lg\t : ChiSqr + a*S\n", ChiSqr + a*S);
			fprintf (path, "%s\t%s\n", "r", "f*dr");
			for (i = 1; i <= N; i++)
				fprintf (path, "%lg\t%lg\n", r[i], x[i]);
			fclose (path);
		}
#endif 
        fflush (stdout);
	}
#endif

	iterationMaximum = (maxIter < ITERATION_MAXIMUM) ? maxIter : ITERATION_MAXIMUM;
	for (i = 1; i <= iterationMaximum; i++) {
		midpoint = (upper + lower) * 0.5;	/* search by bisection */

		Solve (midpoint, I,s,M, r,x,N, G, &ChiSqr, 
                        &S, B,A,D,H,w,d,indx, &a, method);

#ifdef __MWERKS__
		printf ("%15g    %15g    %15g    %15g\n",
			   midpoint,  ChiSqr,    S,    ChiSqr + a*S);
#else
		printf ("%15lg    %15lg    %15lg    %15lg\n",
			   midpoint,  ChiSqr,    S,    ChiSqr + a*S);
#endif 
        fflush (stdout);
		if (ChiSqr > (double) M) {	/* look in lower half */
			upper = midpoint;
		} else {					/* look in upper half */
			lower = midpoint;
		}
		if ( fabs(ChiSqr - (double) M) <= tolerance) break;
	}
	/*
	 * 11. Calculate dx[1..N] from covariance matrix
	 */

	free_ivector (indx, 1, N);
	free_vector (B, 1, N);
	free_vector (w, 1, M);
	free_matrix (A, 1, N, 1, N);
	free_matrix (D, 1, N, 1, N);
	free_matrix (H, 1, N, 1, N);
	return (1 == 1);		/* returns TRUE meaning success */
}


void call_nnls(double **A, int M, int N, double *B, double *X)
{
   double *w;
   int *PZ, result;


//SHOW_VECTOR("B", B, 1, 4);
//SHOW_VECTOR("A[1]", A[1], 1, 4);
//SHOW_VECTOR("A[2]", A[2], 1, 4);
//SHOW_VECTOR("A[3]", A[3], 1, 4);
//SHOW_VECTOR("A[4]", A[4], 1, 4);

	/* Non-Negative Least Squares solution of Ax=B */

    w = vector(1,N);
    PZ = ivector(1,N);
    result = nnls2(A, M, N, B, X, w, PZ);
    free_vector(w, 1,N);
    free_ivector(PZ, 1,N);


//SHOW_VECTOR("X", X, 1, 4);
//SHOW_VECTOR("B", B, 1, 4);

//DEBUG_MARKER;
}
