/* sizes.c
 *
 * General size distribution program for SAS analysis
 */

#ifdef THINK_C
  #include <console.h>
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "recipes.h"
#include "nnls2.h"

#ifdef __MWERKS__
  #include <SIOUX.h>
  #include <console.h>
  /*
   * possibilites:
   *  _dest_os == __win32_os
   *  _dest_os == __mac_os
   *  __MWERKS__
   *  __POWERPC__  _powerc __powerc
   *  macintosh
   */
  #if _dest_os == __win32_os
    #include <WinSIOUX.h>
  #endif
#endif

#define BILLBOARD "sizes (2004.12.23): Analysis of Small-Angle Scattering data"

double  *Q, *SAS, *dSAS;    /* input data */
int     nSAS;               /* number of data points */
double  *r, *dr, *f, *b;    /* solved size distribution */
double  **G;                /* mapping between SAS[] and f[] */

#define MAX_DESC	(4000)	/* allow for up to 4000 character description */
/* the input terms: */
char    ProjectName[256];   /* base name for the output files */
char    SASfile[256];       /* name of the input data file */
char	description[MAX_DESC];	/* for any commented lines in input data */
double  qMin, qMax;         /* units: 1/A */
double  RhoSq;              /* scattering contrast, 10^28, 1/m^4 */
double  fac;                /* 1/cm per data units */
double  err;                /* std. dev. per data units */
double  Bkg;                /* constant to subtract from SAS, data units */
int     shapeModel;			/* 1 = spheres, 2-... not specified yet */
double  Aspect;             /* aspect ratio of scatterer */
double  sLengt;             /* slit length */
int     BinType;            /* linear = 1, logarithmic = 0 */
int     nRadii;             /* number of size distribution points */
double  Dmin, Dmax;         /* units: A */
int     vPower;             /* power of volume weighting in solution */
double  defaultImage;       /* b[], expectation of f[] without constraint */
int     IterMax;            /* maximum number of MaxEnt iterations */
int     iPlot;              /* interval between interim plots */
double  dE_E;               /* dE/E of incident radiation */
int     AnalysisType;       /* 0=reg., 1=MaxEnt, 2=reg+NNLS, 3=NNLS, 4=SVD */
char    timeDate[256];      /* for time/date stamping */
char    *name_template[] =  /* the 3 different types of distribution file */
            {"%s.N-D", "%s.f-D", "%s.i-D"};
char    *column_header[] =  /* column headings for each file type */
            {"N,1/A^4", "f,1/A", "i,A^2"};


#ifndef TRUE
#define TRUE    (1 == 1)
#endif
#ifndef FALSE
#define FALSE   (!TRUE)
#endif

#define REG_LU_METHOD     0         /* standard regularization */
#define MAXENT_METHOD     1         /* Maximum Entropy */
#define REG_NNLS_METHOD   2         /* regularization with NNLS */
#define NNLS_METHOD       3         /* NNLS (no regularization) */
#define SVD_METHOD        4         /* Singular Value Decomposition */
#define SV_CUTOFF         (1.0e-3)

void RecordTheTime (char *s);
int GetInf (char *cmdFile);
int ReadXYZ (char *theFile, int *n, double **x, double **y, double **z);
void MakeSeries (int n, double *x, double *y, 
                double lo, double hi, int lin);
void opus (int N, double *image, int M, double *data);
void tropus (int M, double *data, int N, double *image);
int main (int argc, char **argv);
int Regularize_v2 (double *I, double *s, int M, 
                   double *r, double *x, int N, double **G, 
                   int maxIter, int method);


void RecordTheTime (char *s)
{   /* is this THINK_C specific? */
    time_t now;
    struct tm *date;

    now = time (NULL);
    date = localtime (&now);
        /* example: 92.12.30 22:47:13 */
    strftime (s, 256, "%y.%m.%d, %H:%M:%S", date);
}


#define FGETS_S     fgets (s, 256, path)

int GetInf (char *commandFile)
{
    FILE    *path;
    char    s[256];

    path = fopen (commandFile, "r");
    printf ("Using command file: %s\n", commandFile);
    fflush (stdout);
    if (!path) return FALSE;
    FGETS_S; sscanf (s, "%s", ProjectName);
    FGETS_S; sscanf (s, "%s", SASfile);
    FGETS_S; sscanf (s, "%lg%lg", &qMin, &qMax);
    FGETS_S; sscanf (s, "%lg", &RhoSq);
    FGETS_S; sscanf (s, "%lg", &fac);
    FGETS_S; sscanf (s, "%lg", &err);
    FGETS_S; sscanf (s, "%lg", &Bkg);
    FGETS_S; sscanf (s, "%d", &shapeModel);
    FGETS_S; sscanf (s, "%lg", &Aspect);
    FGETS_S; sscanf (s, "%d", &BinType);
    FGETS_S; sscanf (s, "%d", &nRadii);
    FGETS_S; sscanf (s, "%lg%lg", &Dmin, &Dmax);
    FGETS_S; sscanf (s, "%d", &vPower);
    FGETS_S; sscanf (s, "%lg", &defaultImage);
    FGETS_S; sscanf (s, "%d", &IterMax);
    FGETS_S; sscanf (s, "%lg", &sLengt);
    FGETS_S; sscanf (s, "%lg", &dE_E);
    FGETS_S; sscanf (s, "%d", &AnalysisType);
    fclose (path);
    /* here is an example command file for the test distribution:
                test : Project Name  (only 1st item is read)
            test.sas : SAS file, contains columns: Q  i  esd
    1e-08        100 : qMin qMax, 1/A  (1.0e-8 to 100 means all data)
                 100 : rhosq       : scattering contrast, 10^20 1/cm^-4
                   1 : fac         :   I = fac * ( i - bkg )
                   1 : err         : ESD = fac * err * esd
                 0.1 : bkg         :   I = fac * ( i - bkg )
                   1 : shapeModel  (1=spheroids, no others yet)
                   1 : Aspect Ratio
                   1 : Bin Type    (1=Lin, 0=Log)
                  40 : nRadii
       25        900 : dMin dMax, A
                   1 : n, in N(D)*V^n
              1.0e-6 : defaultDistLevel  (MaxEnt only)
                  32 : IterMax
                   0 : slitLength, 1/A
              0.0002 : dLambda/Lambda
                   0 : analysisType (0=reg., 1=MaxEnt, 2=reg+NNLS, 3=NNLS, 4=SVD)
    */
#ifdef THINK_C
        sprintf (s, "%s.log", ProjectName);
        cecho2file (s, 0, stdout);
#endif
    RecordTheTime (timeDate);
    printf ("\n General Size Distribution "
            "Analysis of SAS data: %s\n\n", timeDate);
    printf ("Input parameters used for this analysis:\n");
#ifdef __MWERKS__
    printf ("%20s : ProjectName\n", ProjectName);
    printf ("%20s : SAS data file\n", SASfile);
    printf ("%9g  %9g : qMin  qMax, 1/A\n", qMin, qMax);
    printf ("%20g : scattering contrast (10^20 cm^-4)\n", RhoSq);
    printf ("%20g : scaling factor\n", fac);
    printf ("%20g : error scaling factor\n", err);
    printf ("%20g : background\n", Bkg);
    printf ("%20d : shape model\n", shapeModel);
    printf ("%20g : aspect ratio, v, as in r * r * v*r\n", Aspect);
    printf ("%20d : bin type (1=Lin, 0=Log)\n", BinType);
    printf ("%20d : number of bins\n", nRadii);
    printf ("%9g  %9g : Dmin & Dmax\n", Dmin, Dmax);
    printf ("%20d : n, in N(D)*V^n\n", vPower);
    printf ("%20g : default distribution level (MaxEnt)\n", defaultImage);
    printf ("%20d : maximum number of MaxEnt iterations\n", IterMax);
    printf ("%20g : slit length, 1/A\n", sLengt);
    printf ("%20g : wavelength distribution, dLambda/Lambda\n", dE_E);
    printf ("%20d : analysisType (0=reg., 1=MaxEnt, 2=reg+NNLS, 3=NNLS, 4=SVD)\n", 
    				AnalysisType);
#else
    printf ("%20s : ProjectName\n", ProjectName);
    printf ("%20s : SAS data file\n", SASfile);
    printf ("%9lg  %9lg : qMin  qMax, 1/A\n", qMin, qMax);
    printf ("%20lg : scattering contrast (10^20 cm^-4)\n", RhoSq);
    printf ("%20lg : scaling factor\n", fac);
    printf ("%20lg : error scaling factor\n", err);
    printf ("%20lg : background\n", Bkg);
    printf ("%20d : shape model\n", shapeModel);
    printf ("%20lg : aspect ratio, v, as in r * r * v*r\n", Aspect);
    printf ("%20d : bin type (1=Lin, 0=Log)\n", BinType);
    printf ("%20d : number of bins\n", nRadii);
    printf ("%9lg  %9lg : Dmin & Dmax\n", Dmin, Dmax);
    printf ("%20d : n, in N(D)*V^n\n", vPower);
    printf ("%20lg : default distribution level (MaxEnt)\n", defaultImage);
    printf ("%20d : maximum number of MaxEnt iterations\n", IterMax);
    printf ("%20lg : slit length, 1/A\n", sLengt);
    printf ("%20lg : wavelength distribution, dLambda/Lambda\n", dE_E);
    printf ("%20d : analysisType (0=reg., 1=MaxEnt, 2=reg+NNLS, 3=NNLS, 4=SVD)\n", 
    				AnalysisType);
#endif 
    printf ("\n");
    fflush (stdout);
    return TRUE;
}


int ReadXYZ (char *theFile, int *n, double **x, double **y, double **z)
{
    FILE    *path;
    int     i, N = 0;
    char    s[256];
    double  *xa, *ya, *za;
    double  xx, yy, zz;     /* for counting lines in the file */
    double  *vector(int nl, int nh);

        /* count the number of lines in the file */
    path = fopen (theFile, "r");
    if (!path) return FALSE;
    while (!feof (path)) {
        fgets (s, 256, path);
        if (s[0] != '#') {
	        i = sscanf (s, "%lg%lg%lg", &xx, &yy, &zz);
	        if (i == 3)  N++;
	    }
    }
    fclose (path);
    /* printf ("There are %d lines in file %s\n", N, theFile);  */
    *n = N;

        /* allocate space (1-based) and read the data */
    xa = *x = vector (1, N);    xa++;
    ya = *y = vector (1, N);    ya++;
    za = *z = vector (1, N);    za++;
    path = fopen (theFile, "r");
    if (!path) return FALSE;

    i = 0;
    description[0] = 0;
    while (!feof (path)) {
        fgets (s, 256, path);
        if (s[0] == '#') {
	    	if (strlen(description) + strlen(s) + 1 < MAX_DESC)
		    	strcat (description, s);
	    } else {
	        i++;
	        sscanf (s, "%lg%lg%lg", xa++, ya++, za++);
	    }
	    if (i == N) break;
    }
    fclose (path);
    return TRUE;
}


void MakeSeries (int n, double *x, double *y, double lo, double hi, int lin)
{
    double mult, step;
    int i;
    if (lin) {              /* algebraic series */
        mult = 1.;
        step = (hi-lo) / (n-1.);
    } else {                /* geometric series */
        mult = exp ( log (hi/lo)/(n-1.) );
        step = 0.;
    }
    x[1] = lo;
    for (i = 2; i <= n; i++)
        x[i] = x[i-1] * mult + step;
    if (y)
        for (i = 1; i <= n; i++)
            y[i] = x[i] * (mult - 1) + step;
}


void opus (int N, double *image, int M, double *data)
/* solution-space -> data-space */
{
    double *sum;
    int i,j;

    for (j = 1; j <= M; j++) {
        sum = data+j;        /* &data[j] */
        for (*sum = 0.0, i = 1; i <= N; i++)
            *sum += image[i] * G[i][j];
    }
}


void tropus (int M, double *data, int N, double *image)
/* data-space -> solution-space */
{
    double *sum;
    int i,j;

    for (i = 1; i <= N; i++) {
        sum = image+i;      /* &image[i] */
        for (*sum = 0.0, j = 1; j <= M; j++)
            *sum += data[j] * G[i][j];
    }
}



int Direct_NNLS_solution (
     double *I, double *s, int M,    /* measured data */
     double *r, double *x, int N,    /* solution vector */
     double **G)                     /* design matrix */
     /*
       NOTE:  M & N definitions are reversed!
       here: M = number of measurements
             N = distribution bins
             G is the transpose of what we need
      */
{
   double *w, **a;
   int *PZ, result, i, j;

   w = vector(1,N);
   PZ = ivector(1,N);
   a = matrix(1,M, 1,N);

   /*
      Since NNLS2 will blow away the design matrix G,
      we copy it to a disposable matrix.
      QUESTION:
        Don't we want to consider the estimated uncertainties?
        Don't we want to consider the bin sizes?
    */
   for (i=1; i<=M; i++) {
     for (j=1; j<=N; j++) {
       a[i][j] = G[j][i];    /* here: a = transpose(G) */
     }
   }
   result = nnls2(a, M, N, I, x, w, PZ);

   free_matrix(a, 1,M, 1,N);
   free_vector(w, 1,N);
   free_ivector(PZ, 1,N);
   return (!result);
}



int Direct_SVD_solution (
     double *I, double *s, int M,    /* measured data */
     double *r, double *x, int N,    /* solution vector */
     double **G,                     /* design matrix */
     double cutoff)                  /* fractional cutoff for singular values */
     /*
       NOTE:  M & N definitions are reversed!
       here: M = number of measurements
             N = distribution bins
             G is the transpose of what we need
      */
{
   double *w, **v, **a, wMax, wMin;
   int *PZ, result, i, j;

   w = vector(1,N);
   v = matrix(1,N, 1,N);
   a = matrix(1,M, 1,N);

   /*
      Since SVD will blow away the design matrix G,
      we copy it to a disposable matrix.
      QUESTION:
        Don't we want to consider the estimated uncertainties?
        Don't we want to consider the bin sizes?
    */
   for (i=1; i<=M; i++) {
     for (j=1; j<=N; j++) {
       a[i][j] = G[j][i];    /* a = transpose(G) */
     }
   }

   /*
       Singular Value decomposition a[][] = u[][] w[] v[][] 
       a[][] is changed into u[][] by SVD
    */
   svdcmp(a, M, N, w, v);
   /* find and zero the insignificant singular values */
   wMax = w[1];
   for (j=1; j<=N; j++) if (w[j]>wMax) wMax = w[j];
   wMin = wMax * cutoff;
   for (j=1; j<=N; j++) if (w[j]<wMin) w[j] = 0.0;

   /*
       Singular Value backsubstitution solution of a z = b
    */
   svbksb(a, w, v, M, N, I, x);

   free_matrix(a, 1,M, 1,N);
   free_matrix(v, 1,N, 1,N);
   free_vector(w, 1,N);
   return (1);
}



int main (int argc, char **argv)
{
    int     proceed, i, j, nIter;
    char    s[256], analysisMethod[256];
    FILE    *path;
    double  *fit;        /* calculation of the SAS fit */
    double  *vector(int nl, int nh);
    void    free_vector (double *v, int nl, int nh);
    double  **matrix(int nrl, int nrh, int ncl, int nch);
    void    free_matrix (double **m, int nrl, int nrh, int ncl, int nch);
    double  Volume (double r, double aspect);
    void    Spheroid (double r[], int nR, double h[], int nH,
                    double **G, double v, double dWave);
    int MaxEnt (
        int n, double *f, double *base, 
        int npt, double *datum, double *sigma,
        double flat, int *iter, int IterMax);


#ifdef __MWERKS__
    SIOUXSettings.asktosaveonclose = 0;
    printf ("\n%s\ninput command file name ==> ", BILLBOARD);
    scanf ("%s", s);
    proceed = GetInf (s);
#else
  #ifdef THINK_C
      argc = ccommand (&argv);
  #endif
    if (argc != 2) {
        printf ("\n%s\nUsage:  sizes  inpFile\n", BILLBOARD);
        printf ("       where inpFile has a special format.\n");
        return 1;
    }
    proceed = GetInf (argv[1]);
#endif 
    printf ("\n%s\n\n", BILLBOARD);
    fflush (stdout);
    if (!proceed) {
        printf ("\n\n Routine GetInf returned FALSE\n");
        return proceed;     /* program ends if GetInf returns FALSE */
    }
    r = vector (1, nRadii);
    dr = vector (1, nRadii);
    f = vector (1, nRadii);
    b = vector (1, nRadii);
    printf ("\n");
    printf ("Reading SAS from file: %s ...", SASfile);
    proceed = ReadXYZ (SASfile, &nSAS, &Q, &SAS, &dSAS);
    if (!proceed) {
        printf ("\n\n Routine ReadXYZ returned FALSE\n");
        return proceed;     /* program ends if ReadXYZ returns FALSE */
    }
    printf ("%d points were read.\n", nSAS);
    fflush (stdout);
    for (i = 1, j=0; i <= nSAS; i++)
        if (qMin <= Q[i] && Q[i] <= qMax) {
            j++;
            Q[j] = Q[i];
            SAS[j] = (SAS[i] - Bkg) * fac;
            dSAS[j] = dSAS[i] * fac * err;
        }
    nSAS = j;
    printf ("%d points remained after Q-range editing.\n", nSAS);
    if (strlen(description) > 0)
    	printf ("file comments:\n%s", description);
    
    if (BinType == 1) {
      printf ("Linear-spaced size bins\n");
    } else {
      printf ("Logarithmic-spaced size bins\n");
    }
    fflush (stdout);
    MakeSeries (nRadii, r, dr, Dmin/2, Dmax/2, BinType);

    RecordTheTime (timeDate);
    printf ("Preparing the G(h,r) matrix ...(%s)\n", timeDate);
    fflush (stdout);
    G = matrix (1, nRadii, 1, nSAS);
    switch (shapeModel) {
    	default:
    	case 1:
            printf ("Spheroidal scatterers ...\n");
    		Spheroid (r, nRadii, Q, nSAS, G, Aspect, dE_E);
    		break;
    }

        /* Do the sizes analysis */
    RecordTheTime (timeDate);
    printf ("Beginning the solution routine ... (%s)\n", timeDate);
    fflush (stdout);
    switch (AnalysisType)  {

        case MAXENT_METHOD: 
             strcpy (analysisMethod, "MaxEnt");
             printf ("%s analysis ...\n", analysisMethod);
             proceed = MaxEnt (nRadii, f, b, nSAS, SAS, dSAS,
                     defaultImage, &nIter, IterMax);
             break;

        case REG_NNLS_METHOD: 
             strcpy (analysisMethod, "regularization with NNLS");
             printf ("%s analysis ...\n", analysisMethod);
             proceed = Regularize_v2 (SAS, dSAS, nSAS,
                     r, f, nRadii, G, IterMax, AnalysisType);
             break;

        case NNLS_METHOD: 
             strcpy (analysisMethod, "NNLS");
             printf ("%s analysis ...\n", analysisMethod);
             proceed = Direct_NNLS_solution (SAS, dSAS, nSAS,
                     r, f, nRadii, G);
             break;

        case SVD_METHOD: 
             strcpy (analysisMethod, "SVD");
             printf ("%s analysis ...\n", analysisMethod);
             proceed = Direct_SVD_solution (SAS, dSAS, nSAS,
                     r, f, nRadii, G, SV_CUTOFF);
             break;

        default:
        case REG_LU_METHOD: 
             strcpy (analysisMethod, "regularization");
             printf ("%s analysis ...\n", analysisMethod);
             proceed = Regularize_v2 (SAS, dSAS, nSAS,
                     r, f, nRadii, G, IterMax, AnalysisType);

    }
    fflush (stdout);

    RecordTheTime (timeDate);
    if (!proceed) {
        printf ("\n\n No solution, %s\n", timeDate);
        if (AnalysisType == 1)
            printf ("No Convergence in %d iterations\n", nIter);
        printf ("Project was: %s\n", ProjectName);
        printf (">>>>>>>>> NO SOLUTION <<<<<<<<<<<\n");
        return proceed;     /* program ends if proceed == FALSE */
    }

        /* write out the size distribution */
    sprintf (s, name_template[vPower], ProjectName);
    path = fopen (s, "w");
    if (path) {
        strcpy (s, column_header[vPower]);
        fprintf (path, 
            "# project: %s, %s analyzed: %s\n",
            ProjectName, analysisMethod, timeDate);
        if (strlen(description) > 0)
	        fprintf (path, "# data file comment:\n%s#\n", description);
        fprintf (path, "%s\t%s\n", "# D,A", s);
        for (i = 1; i <= nRadii; i++)
            fprintf (path, 
#ifdef __MWERKS__
                "%g\t%g\n", 
#else
                "%lg\t%lg\n", 
#endif 
                2*r[i], f[i]/dr[i]/2
            );
        fclose (path);
    }

        /* write out the fitted SAS profile */
    fit = vector (1, nSAS);
    opus (nRadii, f, nSAS, fit);
    sprintf (s, "%s.fit", ProjectName);
    path = fopen (s, "w");
    if (path) {
        fprintf (path, 
            "# project: %s, %s analyzed: %s\n",
            ProjectName, analysisMethod, timeDate);
        if (strlen(description) > 0)
	        fprintf (path, "# data file comment:\n%s#\n", description);
        fprintf (path, 
            "%s\t%s\t%s\t%s\t%s\n",
            "# Q,1/A", "SAS,1/cm", "dSAS,1/cm", "fit,1/cm", "z"
        );
        for (j = 1; j <= nSAS; j++)
            fprintf (path, 
#ifdef __MWERKS__
                "%g\t%g\t%g\t%g\t%f\n", 
#else
                "%lg\t%lg\t%lg\t%lg\t%lf\n", 
#endif 
                Q[j], SAS[j], dSAS[j], fit[j], (SAS[j]-fit[j])/dSAS[j]
            );
    }
    fclose (path);


        /* write out the statistical analysis to stdout */
    statisticalAnalysis( stdout, r,dr,f,nRadii, vPower, Aspect, 
    					Q,SAS,dSAS,fit,nSAS, Bkg, AnalysisType );
    fflush (stdout);

    printf ("\nAnalysis (%s) completed: %s\n", ProjectName, timeDate);

    free_vector (r, 1, nRadii);
    free_vector (dr, 1, nRadii);
    free_vector (f, 1, nRadii);
    free_vector (b, 1, nRadii);
    free_vector (Q, 1, nSAS);
    free_vector (SAS, 1, nSAS);
    free_vector (dSAS, 1, nSAS);
    free_vector (fit, 1, nSAS);
    free_matrix (G, 1, nRadii, 1, nSAS);

    return 0;
}
