/* grid.h
 *
 * Make the grid for the analysis of SAS
 */

/* #define WRITE_INTERMEDIATES      /* useful when debugging */
#define Pi              (3.141592653589793)
#define STEP(a,b,N)     (exp( log(b/a)/( N-1.0 ) ) + 1.0e-5)    /* step-size */
#define SPHERE_APPROX   (0.01)          /* how close to a sphere? */
#define NUM_TABLE_VALS  (500)           /* # of table values */
#define NUM_INTEGRANDS  (100)           /* # of integrands */
#define MIN_DISPERSION  (.005)          /* for wavelength smearing */
#define WaveMean        (1.0)           /* mean wavelength is unity */
#define WaveMin         (WaveMean-dWave)
#define WaveMax         (WaveMean+dWave)
#define TRIANGULAR_DIST (1)
#define GAUSSIAN_DIST   (2)
#define WL_DIST_TYPE    (TRIANGULAR_DIST)

int WriteSizesFile (double *x, double *y, int n, char *fileName);
static double SphereForm (double u);
static double ThisForm (double here, double hor[], double ver[], int n, int j);
double Area (double *x, double *y, int N);
double Volume (double r, double aspect);
double Gaussian(double x, double c, double s);
void Spheroid (double r[], int nR, double h[], int nH,
    double **G, double v, double dWave);
void SpheresVector (double *u, double *S, int nU);
void BinVectorU (double *u, int nU, double v, double uMin, double uMax);
void EllipsoidIntegral (double *u, double *S, int n, 
	double v, double tMin, double tMax);
void WavelengthIntegral (double *u, double *S, int n, 
	double dWave, int wl_dist_type, double tMin, double tMax);
