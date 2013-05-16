/* grid.c
 *
 * Make the grid for the analysis of SAS
 */

#include <stdio.h>
#include <math.h>
#include "grid.h"

extern double RhoSq;
extern int    vPower;
extern double sLengt;

/****************************************
 **************************************** WriteSizesFile
 ****************************************
    Writes the given vectors into a data file 
    (for debugging purposes).
*/
int WriteSizesFile (double *x, double *y, int n, char *fileName)
{
    FILE *path;
    int i;
    path = fopen (fileName, "w");
    if (path) {
        for (i = 1; i <= n; i++)
#ifdef __MWERKS__
            fprintf (path, "%g\t%g\n", x[i], y[i]);
#else
            fprintf (path, "%lg\t%lg\n", x[i], y[i]);
#endif 
        fclose (path);
    }
    return (0);
}

/****************************************
 **************************************** SphereForm
 ****************************************
    Calculates the form factor for scattering from
    a single, homogeneous sphere according to Debye.
*/
static double SphereForm (double u)
{               /* standard form factor for spheres */
    double F;
    F = 3.0 * (sin(u) - u*cos(u)) / (u*u*u);
    return F*F;     /* square it for intensity */
}


static double ThisShape (double here, double hor[], double ver[], 
		int n, int j)
{
    double   x1, y1,   x2, y2, result;
    if (j > n)
      return SphereForm (here); /*  no table data so get exactly */
    if (j == n)
        j--;                    /* make an extrapolation */

    /* otherwise, use 2-pt. linear interpolation */
    x1 = hor[j];
    y1 = ver[j];
    x2 = hor[j+1];
    y2 = ver[j+1];
    result = y1 * (here - x1) / (x2 - x1)
           + y2 * (x2 - here) / (x2 - x1);

#ifdef WRITE_INTERMEDIATES
    if (result <=0.0)
#ifdef __MWERKS__
        printf("non-positive shape function!  S(%g)=%g, j=%d, n=%d\n",
                here, result, j, n);
#else
        printf("non-positive shape function!  S(%lg)=%lg, j=%d, n=%d\n",
                here, result, j, n);
#endif 
#endif
    return result;
}

double Area (double *x, double *y, int N)
{                       /* integration by trapezoid rule */
    double sum = 0.0;
    int i;
    for (i = 1; i < N; i++)
        sum += (x[i+1] - x[i]) * (y[i+1] + y[i]);
    return sum*0.5;
}

#define Root2Pi (2.506628274631000502)
double Gaussian(double x, double c, double s)
{
    double g = (x - c) / s;
    return (1 / s / Root2Pi) * exp(-0.5*g*g);
}


double Volume (double r, double aspect)
{
    double vol = Pi * 4 / 3 * aspect * pow((double) r*1.e-8, (double) 3);
    /*  printf ("Volume: r=%g aspect=%g volume=%g\n", r, aspect, vol); */
    return vol;
}


/****************************************
 **************************************** SpheresVector
 ****************************************
 *
 * Calculate the form factor for spheres
 *
 */
void SpheresVector (double *u, double *S, int n)
{
	int i;
	for (i = 1; i <= n; i++) {
		S[i] = SphereForm (u[i]);
		/* printf ("SpheresVector: %d/%d  u=%g  S=%g\n", i, n, u[i], S[i]); */
	}
	return;
}


/****************************************
 **************************************** BinVectorU
 ****************************************
 *
 * Bin the vector u[*]
 *
 */
void BinVectorU (double *u, int n, double v, double uMin, double uMax)
{
	int i;
	double vMin, vMax, du;

    vMin = (v < 1.0) ? v : 1.0;
    vMax = (v < 1.0) ? 1.0 : v;
    uMin *=  0.95 * vMin;		/* expand the limits just a skosh */
    uMax *=  1.05 * vMax;
    du = STEP(uMin, uMax, n);   /* the step size */
    u[1] = uMin;
    for (i = 2; i <= n; i++)
      u[i] = u[i-1] * du;

	return;
}


/****************************************
 **************************************** EllipsoidIntegral
 ****************************************
 *
 * Integrate to get the ellipsoidal form factor
 *
 */
void EllipsoidIntegral (
	double *u, 
	double *S, 
	int n, 
	double v, 
	double tMin, 
	double tMax)
{
	int i, j, index;
	double *t, *SS, *x, *y, dt, uNow, tNow;
    void hunt(double xx[], int n, double x, int *jlo);
    double  *vector(int nl, int nh);
    void    free_vector (double *v, int nl, int nh);

    t = vector (1, n);
    SS = vector (1, n);
    x = vector (1, NUM_INTEGRANDS);
    y = vector (1, NUM_INTEGRANDS);
    /* 
     * ordinate of integration 
     * (expand the limits just a skosh)
     */
    tMin *= 0.95;
    tMax *= 1.05;
    for (j = 1; j <= NUM_INTEGRANDS; j++)
    	x[j] = (double) (j - 1.0) / (NUM_INTEGRANDS - 1.0);
    dt = STEP(tMin, tMax, n);						/* the step size */
    tNow = tMin;                                    /* first point */
    /*
     * for each point in the vector S[*] ...
     */
    for (i = 1; i <= n; i++) {
    	t[i] = tNow;                                /* keep for lookup table */
    	/*
    	 * ... set up the integrand for numerical integration at tNow ...
    	 */
    	for (j = 1; j <= NUM_INTEGRANDS; j++) {
          uNow = tNow * sqrt (1 + x[j]*x[j] * (v*v - 1));
          hunt (u, n, uNow, &index);                /* get the index */
          if (index == 0) index = 1;                /* just in case */
          y[j] = ThisShape(uNow, u,S,n, index);     /* get the value */
    	}
    	/*
    	 * ... and integrate
    	 */
    	SS[i] = Area(x, y, NUM_INTEGRANDS);         /* integration */
    	tNow = tNow * dt;                           /* for the next t */
	}
	/*
	 * copy the results back into the input vectors
	 */
	for (i = 1; i <= n; i++) {
		u[i] = t[i];	/* new ordinates! */
		S[i] = SS[i];
	}

    free_vector (t, 1, n);
    free_vector (SS, 1, n);
    free_vector (x, 1, NUM_INTEGRANDS);
    free_vector (y, 1, NUM_INTEGRANDS);
	return;
}


/****************************************
 **************************************** WavelengthIntegral
 ****************************************
 *
 * Integrate to include the wavelength distribution
 *
 */
void WavelengthIntegral (
	double *u, 
	double *S, 
	int n, 
	double dWave,
	int wl_dist_type,
	double tMin, double tMax)
{
	int i, j, index=1;
	char msg[256];
	double *p, *x, *y, dw, *t, *SS, dt, tNow, uNow, WaveNorm;
    void hunt(double xx[], int n, double x, int *jlo);
    double  *vector(int nl, int nh);
    void    free_vector (double *v, int nl, int nh);

	/*
	 * define the wavelength distribution, p(x)
	 */
    p = vector (1, NUM_INTEGRANDS);
    x = vector (1, NUM_INTEGRANDS);
	dw = (double) (WaveMax - WaveMin) / (NUM_INTEGRANDS - 1.0);
	for (j = 1; j <= NUM_INTEGRANDS; j++) {
		/*
		 * wavelengths for integration
		 */
		x[j] = WaveMin + (j - 1.0) * dw;
		switch (wl_dist_type) {
			default:
			case 1:					/* triangular distribution function */
				p[j] = ((double) 1.0 - fabs((x[j]-WaveMean)/dWave))/dWave;
				break;

			case 2:					/* Gaussian distribution function */
				p[j] = Gaussian (x[j], WaveMean, dWave);
		}
	}
	WaveNorm = Area(x, p, NUM_INTEGRANDS);        /* to make unit area */
#ifdef WRITE_INTERMEDIATES
#ifdef __MWERKS__
	printf (" dWave = %g\n WaveNorm = %g\n Wave step = %g\n", 
                dWave, WaveNorm, dw);
#else
	printf (" dWave = %lg\n WaveNorm = %lg\n Wave step = %lg\n", 
                dWave, WaveNorm, dw);
#endif 
	WriteSizesFile (x, p, NUM_INTEGRANDS, "P_lambda.dat");
#endif

    t = vector (1, n);
    SS = vector (1, n);
    y = vector (1, NUM_INTEGRANDS);
    dt = STEP(tMin, tMax, n);    /* the step size */
    tNow = tMin;                  /* first point */
    /*
     * For every data point, t[i] ...
     */
    for (i = 1; i <= n; i++) {
    	t[i] = tNow;
    	/*
    	 * ... evaluate the integrand ...
    	 */
    	for (j = 1; j <= NUM_INTEGRANDS; j++) {
    		uNow = tNow * WaveMean / x[j];
    		hunt (u, n, uNow, &index);
            if (index == 0) index = 1;                /* just in case */
    		y[j] = p[j] * ThisShape(uNow, u,S,n, index);
    	}
#ifdef WRITE_INTERMEDIATES
    	printf ("Writing file: wave-%03d.dat\n", i);
    	if (!(i % 5000)) {
    		sprintf (msg, "wave-%03d.dat", i);
    		WriteSizesFile (x, y, NUM_INTEGRANDS, msg);
    	}
#endif
		/*
		 * ... end integrate
		 */
    	SS[i] = Area(x, y, NUM_INTEGRANDS) / WaveNorm;
    	tNow *= dt;               /* for the next t */
    }

	/*
	 * copy the results back into the input vectors
	 */
	for (i = 1; i <= n; i++) {
		u[i] = t[i];	/* new ordinates! */
		S[i] = SS[i];
	}

    free_vector (t, 1, n);
    free_vector (SS, 1, n);
    free_vector (p, 1, NUM_INTEGRANDS);
    free_vector (x, 1, NUM_INTEGRANDS);
    free_vector (y, 1, NUM_INTEGRANDS);
	return;
}


/****************************************
 **************************************** SlitSmear
 ****************************************
 *
 * Integrate to include slit smearing
 *
 */
double SlitSmear (
	double *u, 
	double *S,
	int n,
	double slitLength,
	double h,
	double r)
{
	int i, index = 1;
	double *x, *y, dx, xNow, uNow, shape, p;
    void hunt(double xx[], int n, double x, int *jlo);
    double  *vector(int nl, int nh);
    void    free_vector (double *v, int nl, int nh);

    x = vector (1, NUM_INTEGRANDS);
    y = vector (1, NUM_INTEGRANDS);
    /* 
     * set the step size and first point
     */
    dx = slitLength / (double) (NUM_INTEGRANDS - 1);
    xNow = 0.0;
    for (i = 1; i <= NUM_INTEGRANDS; i++) {
    	x[i] = xNow;
    	p = (xNow < slitLength) ? 1.0 : 0.0;	/* rectangular slit */
    	uNow = r * sqrt( h*h + xNow*xNow );		/* circularly-symmetric */
    	hunt (u, n, uNow, &index);
        if (index == 0) index = 1;              /* just in case */
    	y[i] = p * ThisShape(uNow, u,S,n, index);
    	xNow += dx;
    }
    shape = (1/slitLength) * Area(x, y, NUM_INTEGRANDS);
    free_vector (x, 1, NUM_INTEGRANDS);
    free_vector (y, 1, NUM_INTEGRANDS);
    return shape;
}


/****************************************
 **************************************** Spheroid
 ****************************************

    Particle form factor for spheroids of R x R x vR, R = radius
      ref: Roess and Shull; J Appl Phys 18 (1947) 308-313.
      coded: PRJ on 7.2.90
      modified: PRJ on 4.May.92 included triangular wavelength dispersion
    The method of eq.6 (second part) is implemented by direct evaluation
      using eq. (4) & (5) and numerical averaging.
    Lookup tables are used to speed up the calculations.

*/
void Spheroid (
    double r[], int nR,  /* array of scattering vectors numbered 1 to nH */
    double h[], int nH,  /* array of radii numbered 1 to nR */
    double **G,  /* grid of (form factors)**2 numbered 1 to nH & 1 to nR */
    double v,            /* aspect ratio */
    double dWave)        /* FWHM of wavelength */
{
    int     i, j, index;
    char    msg[256];
    double  *u, *S, *SS, weight, hrNow, shape;
    void    hunt(double xx[], int n, double x, int *jlo);
    double  *vector(int nl, int nh);
    void    free_vector (double *v, int nl, int nh);
    double  uMin, uMax;

    u = vector (1, NUM_TABLE_VALS);
    S = vector (1, NUM_TABLE_VALS);
    SS = vector (1, NUM_TABLE_VALS);

	/*
	 * Calculate the table of "spheres" form factors: S(u)
	 */
	uMin = h[1] * r[1] * (WaveMean/WaveMax);
	uMax = sqrt(h[nH]*h[nH] + sLengt*sLengt) * r[nR] * (WaveMean/WaveMin);
	BinVectorU (u, NUM_TABLE_VALS, v, uMin, uMax );
    SpheresVector (u, S, NUM_TABLE_VALS);
#ifdef WRITE_INTERMEDIATES
    printf ("Writing file: S_s.dat\n");
    WriteSizesFile (u, S, NUM_TABLE_VALS, "S_s.dat");
#endif

/*
 * Integrate to get the ellipsoidal form factor
 * (only integrate if not spherical)
 */
    if (fabs(v-1) > SPHERE_APPROX) {
        printf ("Ellipsoidal scatterers ...\n");
		uMin = h[1] * r[1] * (WaveMean/WaveMax);
		uMax = sqrt(h[nH]*h[nH] + sLengt*sLengt) * r[nR] * (WaveMean/WaveMin);
		EllipsoidIntegral (u, S, NUM_TABLE_VALS, v, uMin, uMax );
#ifdef WRITE_INTERMEDIATES
    	printf ("Writing file: S_e.dat\n");
        WriteSizesFile (u, S, NUM_TABLE_VALS, "S_e.dat");
#endif
    fflush (stdout);
    }


/*
 *  Integrate to include the wavelength distribution
 * (only integrate if finite dispersion)
 */
    if (dWave > MIN_DISPERSION) {
        printf ("Wavelength smearing effects ...\n");
		uMin = h[1] * r[1];
		uMax = sqrt(h[nH]*h[nH] + sLengt*sLengt) * r[nR];
    	WavelengthIntegral (u,S,NUM_TABLE_VALS, dWave, WL_DIST_TYPE, 
    						uMin, uMax );
#ifdef WRITE_INTERMEDIATES
    	printf ("Writing file: S_w.dat\n");
        WriteSizesFile (u, S, NUM_INTEGRANDS, "S_w.dat");
#endif
    }


    if (sLengt < h[1]) {
      printf ("Pinhole collimation ...\n");
    } else {
      printf ("Slit length collimation ...\n");
    }
    fflush (stdout);
	/*
	 *  Calculate G( h[j], r[i] ) from the table
	 */
    for  (i = 1; i <= nR; i++) {        /* at each scattering vector */
      weight = pow ((double) Volume (r[i], v), (double) 2-vPower);
      for  (j = 1; j <= nH; j++) {      /* at each radius */
        hrNow = r[i] * h[j];                /* indexed into table */
		/*
		 *  Include slit smearing effects, if significant
		 */
    	if (sLengt < h[1]) {
    		/*
    		 * perfect collimation
    		 *
    		 * find <shape> by interpolation in the table of S(u) 
    		 */
	        hunt (u, NUM_TABLE_VALS, hrNow, &index);
            if (index == 0) index = 1;                /* just in case */
	        shape = ThisShape(hrNow, u,S,NUM_TABLE_VALS, index);
    	} else {						/* assume slit-smeared */
    		/*
    		 * assume slit-smeared collimation
    		 *
    		 * find <shape> by integration of S(r*sqrt[h^2+l^2]) 
    		 */
    		shape = SlitSmear (u,S,NUM_TABLE_VALS, sLengt, h[j], r[i]);
    	}
        G[i][j] = RhoSq*1.e20 * weight * shape;
      }
    }

    free_vector (u, 1, NUM_TABLE_VALS);
    free_vector (S, 1, NUM_TABLE_VALS);
    free_vector (SS, 1, NUM_TABLE_VALS);
}
