/* maxent.c
 *
 * Entropy maximization routine as described in the article
 * J Skilling and RK Bryan; Mon Not R Astr Soc 211 (1984) 111 - 124.
 * ("Mon Not R Astr Soc" is the "Monthly Notices of the 
 *      Royal Astronomical Society")
 *
 * This program is a translation from the FORTRAN version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TEST_LIMIT      ( 0.05 )    /* for convergence */
#define CHI_SQR_LIMIT   ( 0.01 )    /* maximum difference in ChiSqr for a solution */
#define SEARCH_DIRECTIONS   (3)     /* <10.  Best not to change it now */
#define SD1         (SEARCH_DIRECTIONS+1)
#define RESET_STRAYS    ( 1 )   /* was 0.001, correction of stray negative values */

    double chisq, chtarg, chizer, fSum, blank;
    double beta[SD1], c1[SD1], c2[SD1][SD1], s1[SD1], s2[SD1][SD1];
    double *ox, *z, *cgrad, *sgrad, **xi, **eta;
    static char msg[256];                   /* generic messages */

int MaxEnt (
    int n, double *f, double *base, 
    int npt, double *datum, double *sigma,
    double flat, int *iter, int IterMax);
static void MaxEntMove (int m);
static double Dist (int m);
static double ChiNow (double ax, int m);
static void ChoSol (double a[SD1][SD1], double b[SD1], int n, double *beta);
static void errorMaxEnt (char *error_text);

#ifndef TRUE
#define TRUE    (1 == 1)
#endif
#ifndef FALSE
#define FALSE   (!TRUE)
#endif


#define DITCH_VECTORS   \
    free_vector (ox,    1, npt);                        \
    free_vector (z,     1, npt);                        \
    free_vector (cgrad, 1, n);                          \
    free_vector (sgrad, 1, n);                          \
    free_matrix (xi,    1, SEARCH_DIRECTIONS, 1, n);    \
    free_matrix (eta,   1, SEARCH_DIRECTIONS, 1, npt);  \

int MaxEnt (
    int n, double *f, double *base, 
    int npt, double *datum, double *sigma,
    double flat, int *iter, int IterMax)
{
    int i, j, k, l;
    double exp1, a, b, test, snorm, cnorm, tnorm;
    double S, fChange, df;
    void opus (int N, double *f, int M, double *F);
    void tropus (int M, double *F, int N, double *f);
    double  *vector(int nl, int nh);
    void    free_vector (double *v, int nl, int nh);
    double  **matrix(int nrl, int nrh, int ncl, int nch);
    void    free_matrix (double **m, int nrl, int nrh, int ncl, int nch);

    ox      = vector(1, npt);
    z       = vector(1, npt);
    cgrad   = vector(1, n);
    sgrad   = vector(1, n);
        /* 
         * Note that the order of subscripts for
         * "xi" and "eta" has been reversed from
         * the convention used in the FORTRAN version
         * to enable parts of them to be passed as
         * as vectors to "opus" and "tropus".
         */
    xi      = matrix(1, SEARCH_DIRECTIONS, 1, n);
    eta     = matrix(1, SEARCH_DIRECTIONS, 1, npt);

    exp1 = exp (1.0);
    blank = flat;       /* copy into a local "COMMON" variable */
        /* determine the reference line */
    if (blank == 0.) {
        for (i = 1; i <= n; i++)
            blank += base[i];
        blank /= n;
#ifdef __MWERKS__
        printf ("%g = Average of BASE\n", blank);
#else
        printf ("%lg = Average of BASE\n", blank);
#endif 
    } else {
#ifdef __MWERKS__
        printf ("%g : Constant BASE setting\n", blank);
#else
        printf ("%lg : Constant BASE setting\n", blank);
#endif 
        for (i = 1; i <= n; i++)
            base[i] = blank;
    }

    printf ("\nMaxEnt routine beginning ...\n");
    printf ("%d data points, %d image points\n", npt, n);
        printf ("%3s/%3s", "itr", "iMx");
        printf (" %5s  %8s", "align", "entropy");
        printf (" %12s %10s", "sqrt(X2/N)", "off");
        printf ( "%12s %8s\n", "fSum", "change");
    fflush (stdout);

    chizer = npt;
    chtarg = (chizer = npt);
    for (i = 1; i <= n; i++)
        f[i] = base[i];     /* initial distribution is featureless */

    for (*iter = 1; *iter <= IterMax; (*iter)++) {
        opus (n, f, npt, ox);   /* calc. the model intensity from "f" */
        for (j = 1, chisq = 0.0; j <= npt; j++) {
            a = ox[j] - datum[j];
            a /= sigma[j];
            chisq += a*a;
            ox[j] = 2. * a / sigma[j];
        }
        tropus (npt, ox, n, cgrad);     /* cgrad[i] = del(C)/del(f[i]) */
        test = 0.0;     /* mismatch between entropy and ChiSquared gradients */
        snorm = 0.0;    /* entropy term */
        cnorm = 0.0;    /* ChiSqr term  */
        tnorm = 0.0;    /* norm for the gradient term TEST  */
        fSum = 0.0;     /* find the sum of the f-vector */
        for (i = 1; i <= n; i++) {
            fSum += f[i];
                        /* sgrad[i] = del(S)/del(f[i]) */
            sgrad[i] = -log(f[i]/base[i]) / (blank*exp1);
            snorm += f[i] * sgrad[i]*sgrad[i];      /* Eq. 22 */
            cnorm += f[i] * cgrad[i]*cgrad[i];
            tnorm += f[i] * sgrad[i] * cgrad[i];
        }
        snorm = sqrt(snorm);
        cnorm = sqrt(cnorm);
        a = 1.0;
        b = 1.0 / cnorm;
        if (*iter > 1) {
          test = sqrt( ( 1.0 - tnorm/(snorm*cnorm) )/2 );
          a = 0.5 / (snorm * test);
          b *= 0.5 / test;
        }
        for (i = 1; i <= n; i++) {
            xi[1][i] = f[i] * cgrad[i] / cnorm;
            xi[2][i] = f[i] * (a * sgrad[i] - b * cgrad[i]);
        }
        opus (n, xi[1], npt, eta[1]);   /* image --> data */
        opus (n, xi[2], npt, eta[2]);   /* image --> data */
        for (j = 1; j <= npt; j++)
            ox[j] = eta[2][j] / (sigma[j]*sigma[j]);
        tropus (npt, ox, n, xi[3]);     /* data --> image */
        a = 0.0;
        for (i = 1; i <= n; i++) {
          b = f[i] * xi[3][i];
          a += b * xi[3][i];
          xi[3][i] = b;
        }
        a = 1.0 / sqrt (a);
        for (i = 1; i <= n; i++)
            xi[3][i] *= a;
        opus (n, xi[3], npt, eta[3]);   /* image --> data */

        for (k = 1; k <= SEARCH_DIRECTIONS; k++) {
            s1[k] = 0.0;
            c1[k] = 0.0;
            for (i = 1; i <= n; i++) {
                s1[k] += xi[k][i] * sgrad[i];
                c1[k] += xi[k][i] * cgrad[i];
            }
            c1[k] /= chisq;
        }
        for (k = 1; k <= SEARCH_DIRECTIONS; k++) {
            for (l = 1; l <= k; l++) {
                s2[k][l] = 0.0;
                c2[k][l] = 0.0;
                for (i = 1; i <= n; i++)
                    s2[k][l] -= xi[k][i] * xi[l][i] / f[i];
                for (j = 1; j <= npt; j++)
                    c2[k][l] += eta[k][j] * eta[l][j]
                                / (sigma[j]*sigma[j]);
                s2[k][l] /= blank;
                c2[k][l] *= 2. / chisq;
            }
        }
            /* reflect across the body diagnonal */
        c2[1][2] = c2[2][1];
            c2[1][3] = c2[3][1];
                c2[2][3] = c2[3][2];
        s2[1][2] = s2[2][1];
            s2[1][3] = s2[3][1];
                s2[2][3] = s2[3][2];
    
        beta[1] = -0.5 * c1[1] / c2[1][1];
        beta[2] = 0.0;
        beta[3] = 0.0;
        if (*iter > 1) MaxEntMove(SEARCH_DIRECTIONS);
    
    /*  Modify the current distribution (f-vector)  */
        fSum = 0.0;     /* find the sum of the f-vector */
        fChange = 0.0;      /* and how much did it change?  */
        for (i = 1; i <= n; i++) {
            df = 0.0;
            for (j = 1; j <= SEARCH_DIRECTIONS; j++)
                df += beta[j] * xi[j][i];
            /*
             * As mentioned at the top of p.119,
             * need to protect against stray negative values.
             * In this case, set them to RESET_STRAYS * base[i]
             */
            if ( f[i]+df <= 0.0 )
                df = RESET_STRAYS * base[i] - f[i];
          f[i] += df;                   /* adjust the f-vector  */
          fChange += df;
          fSum += f[i];
        }

            /* calculate the normalized entropy */
        S = 0.0;
        for (i = 1; i <= n; i++)
            S -= (f[i] / fSum) * log (f[i] / fSum);     /* from Skilling and Bryan, eq. 1 */
        opus (n, f, npt, z);            /* image --> data */
        chisq = 0.0;                    /* get the new ChiSquared   */
        for (j = 1; j <= npt; j++) {
            z[j] = (datum[j] - z[j])/sigma[j];  /* standardized residuals */
            chisq += z[j]*z[j];                     /* report this ChiSq */
        }

/*
 *       printf ("\n%d of %d iterations, %d data points, %d image points\n",
 *           *iter, IterMax, npt, n);
 *       printf ("%5.2lf%% misalignment, normalized S = %lg\n", 
 *           100*test, S);
 *       printf (
 *           "%17s:%8s = %12.5lg%10s = %10.4lf\n", 
 *           "sqrt((Chi^2)/n)",
 *           "target",   sqrt (chtarg/npt),
 *           "% off",    (*iter > 1) ? 100*( sqrt(chisq/chtarg)-1) : 0
 *       );
 *       printf (
 *           "%17s:%8s = %12.6lg%10s = %10.4lf\n", 
 *           "f-vector",
 *           "sum",      fSum,
 *           "% change", 100*fChange/fSum
 *       );
 */
#ifdef __MWERKS__
        printf ("%3d/%3d", *iter, IterMax);
        printf (" %5.2f%% %8g", 100*test, S);
        printf (" %12.5g %10.4f", sqrt (chtarg/npt),
            (*iter > 1) ? 100*( sqrt(chisq/chtarg)-1) : 0
        );
        printf ( "%12.6g %8.2f\n", fSum, 100*fChange/fSum);
#else
        printf ("%3d/%3d", *iter, IterMax);
        printf (" %5.2lf%% %8lg", 100*test, S);
        printf (" %12.5lg %10.4lf", sqrt (chtarg/npt),
            (*iter > 1) ? 100*( sqrt(chisq/chtarg)-1) : 0
        );
        printf ( "%12.6lg %8.2lf\n", fSum, 100*fChange/fSum);
#endif 
        fflush (stdout);

            /* See if we have finished our task. */
            /* do the hardest test first */
        if (fabs(chisq/chizer-1.0) < CHI_SQR_LIMIT)
            if (test < TEST_LIMIT) {
                DITCH_VECTORS;
                return TRUE;        /* solution FOUND returns here */
            }
    }

    DITCH_VECTORS;
    return FALSE;       /* solution NOT FOUND returns here */
}


#define MAX_MOVE_LOOPS  (500)     /* for no solution in routine: move */
#define MOVE_PASSES (0.001)  /* convergence test in routine: move */

static void MaxEntMove (int m)
/* "m" is not really a variable.  It is "SEARCH_DIRECTIONS" */
{
    double a1, a2, anew, cmin, ctarg, f1, f2, fx, w;
    int loop, k;

    a1 = 0;         /* lower bracket  "a"   */
    a2 = 1;         /* upper bracket of "a" */
    cmin = ChiNow (a1, m);
    ctarg = (cmin*chisq > chizer) 
        ? (1.0 + cmin)/2 
        : chizer/chisq;
    f1 = cmin - ctarg;
    f2 = ChiNow (a2,m) - ctarg;
    for (loop = MAX_MOVE_LOOPS; loop; loop--) {
        anew = (a1+a2)/2;       /* choose a new "a" */
        fx = ChiNow (anew, m) - ctarg;
        if (f1*fx > 0) {
            a1 = anew;
            f1 = fx;
        }
        if (f2*fx > 0) {
            a2 = anew;
            f2 = fx;
        }
        if (fabs(fx) < MOVE_PASSES) break;
    }

    /*
     *  If the preceding loop finishes,
     *    then we do not seem to be converging.
     *  Stop gracefully because not every computer
     *    uses control-C (etc.) as an exit procedure.
     */
    if (!loop) {
        sprintf (msg, "%s = %d\n%s\n\n%s\n",
            "Loop counter", MAX_MOVE_LOOPS,
            "No convergence in alpha chop (routine: MaxEntMove).",
            "Program cannot continue.");
        errorMaxEnt (msg);
    }

    w = Dist (m);
    if (w > 0.1*fSum/blank) /* invoke the distance penalty? */
        for (k = 1; k <= m; k++)
            beta[k] *= sqrt (fSum/(blank*w));
    chtarg = ctarg * chisq;
}


static double Dist (int m)
/* "m" is not really a variable.  It is "SEARCH_DIRECTIONS" */
{
    double w, z;
    int k, l;
    
    for (k = 1, w = 0; k <= m; k++) {
        for (l = 1, z = 0; l <= m; l++)
            z -= s2[k][l] * beta[l];
        w += beta[k] * z;
    }
    return w;
}


static double ChiNow (double ax, int m)
/* "m" is not really a variable.  It is "SEARCH_DIRECTIONS" */
{
    int k, l;
    double bx = 1 - ax, w, z;
    double a[SD1][SD1], b[SD1];

    for (k = 1; k <= m; k++) {
        for (l = 1; l <= m; l++)
            a[k][l] = bx * c2[k][l] - ax * s2[k][l];
        b[k] = -( bx * c1[k] - ax * s1[k] );
    }
    ChoSol (a, b, m, beta);
    for (k = 1, w = 0; k <= m; k++) {
        for (l = 1, z = 0; l <= m; l++)
            z += c2[k][l] * beta[l];
        w += beta[k] * (c1[k] + z/2);
    }
    return 1 +  w;
}


static void ChoSol (double a[SD1][SD1], double b[SD1], int n, double *beta)
{
    int i, j, k;
    double z, fl[SD1][SD1], bl[SD1];

    if (a[1][1] <= 0) {
#ifdef __MWERKS__
        sprintf (msg, "%s\n%s = %g\n\n%s\n",
            "Fatal error in routine ChoSol.",
            "a[1][1]", a[1][1],
            "Program cannot continue.");
#else
        sprintf (msg, "%s\n%s = %lg\n\n%s\n",
            "Fatal error in routine ChoSol.",
            "a[1][1]", a[1][1],
            "Program cannot continue.");
#endif 
        errorMaxEnt (msg);
    }
    fl[1][1] = sqrt (a[1][1]);
    for (i = 2; i <= n; i++) {
        fl[i][1] = a[i][1] / fl[1][1];
        for (j = 2; j <= i; j++) {
            for (k = 1, z = 0;  k < j; k++)
                z += fl[i][k] * fl[j][k];
            z = a[i][j] - z;
            fl[i][j] = (j == i) ? sqrt(z) : z / fl[j][j];
        }
    }

    bl[1] = b[1] / fl[1][1];
    for (i=2; i <= n; i++) {
        for (k = 1, z = 0;  k < i; k++)
            z += fl[i][k] * bl[k];
        bl[i] = (b[i] - z) / fl[i][i];
    }
    beta[n] = bl[n] / fl[n][n];
    for (i = n - 1; i >= 1; i--) {
        for (k = i+1, z = 0; k <= n; k++)
            z += fl[k][i] * beta[k];
        beta[i] = (bl[i] - z) / fl[i][i]; 
    }
}

static void errorMaxEnt (char *error_text) /* re: NR standard error handler */
{
    fprintf (stderr, "\n");
    fprintf (stderr, "MaxEnt package run-time error...\n%s\n", error_text);
    fprintf (stderr, "... now exiting to system ...\n");
    /*  WaitMouseClick ();  */
    exit (1);   /* return a non-zero value */
}
