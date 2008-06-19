/* stats.c
 *
 * Statistical analysis of size distribution and data fit.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Pi              (3.141592653589793)

int statisticalAnalysis (
    FILE *path, 
    double *r,
    double *dr,
    double *dist,
    int nRadii, 
    int vPower,
    double Aspect,
    double *Q,
    double *SAS,
    double *dSAS,
    double *fit,
    int nSAS,
    double background,
    int AnalysisType);

static void distributionAnalysis (
	double *r, 
	double *f, 
	double *dr, 
	int N, 
	int rPower,
	double Aspect,
	double *area, 
	double *mode, 
	double *mean, 
	double *stdDev,
	int AnalysisType,
	double *entropy);


#ifndef TRUE
#define TRUE    (1 == 1)
#endif
#ifndef FALSE
#define FALSE   (!TRUE)
#endif

	/*
	 * Skilling & Bryan, eq. 1
	 */
#define ENTROPY(FRAC) (- (FRAC) * log (FRAC))

static void distributionAnalysis (double *r, 
	double *f, double *dr, int N, int rPower, double Aspect,
	double *area, double *mode, double *mean, double *stdDev,
	int AnalysisType, double *entropy)
{
	int i, modeIndex;
	double dist, R, rVol, modeVal = 0;
	double sumY = 0.0, sumXY = 0.0, sumXXY = 0.0;
	*area = 0.0;
	for (i = 1; i <= N; i++) {
			/*
			 * calculate the volume for an ellipsoid or radius r[i]
			 */
		R = r[i] * 1.e-8;	/* convert to cm */
		rVol = (4 * Pi / 3) * R*R*R*Aspect;
		dist = f[i] / dr[i] * pow( (double) rVol, (double) rPower);
		if ( (i == 1) || (dist >= modeVal) ) {
			modeIndex = i;
			modeVal = dist;
		}
		*area += dist * dr[i];
		sumY += dist * dr[i]; 
		sumXY += dist * dr[i] * r[i]; 
		sumXXY += dist * dr[i] * r[i] * r[i];
	}
	*mode = r[modeIndex];
	*mean = sumXY / sumY;
	*stdDev = sqrt( (sumXXY / sumY) - pow((double) sumXY / sumY, (double) 2) );
	*entropy = 0.0;
	if (AnalysisType == 1)
	  for (i = 1; i <= N; i++) {
			/*
			 * calculate the volume for an ellipsoid of radius r[i]
			 */
		R = r[i] * 1.e-8;	/* convert to cm */
		rVol = (4 * Pi / 3) * R*R*R*Aspect;
		dist = f[i] / dr[i] * pow( (double) rVol, (double) rPower);
		*entropy += ENTROPY(dist * dr[i] / (*area));
	  }
}

int statisticalAnalysis (
    FILE *path, 
    double *r, double *dr, double *f, int N, 
    int vPower, double Aspect,
    double *Q, double *SAS, double *dSAS, double *fit, int M,
    double background, int AnalysisType)
{
	int    i, j;
	double	nMode,		vMode,		iMode,
			nArea,		vArea,		iArea,
			nMean, 		vMean,		iMean,
			nStdDev, 	vStdDev, 	iStdDev,
			nEntropy, 	vEntropy, 	iEntropy,
			suggested_background, weight, sum1, sum2, background_esd,
			Qsqr_dQ, invariant, invariant_esd, invariant_fit;

	distributionAnalysis (r, f, dr, N, 0 - vPower, Aspect,
		&nArea, &nMode, &nMean, &nStdDev, AnalysisType, &nEntropy);
	distributionAnalysis (r, f, dr, N, 1 - vPower, Aspect,
		&vArea, &vMode, &vMean, &vStdDev, AnalysisType, &vEntropy);
	distributionAnalysis (r, f, dr, N, 2 - vPower, Aspect,
		&iArea, &iMode, &iMean, &iStdDev, AnalysisType, &iEntropy);

	sum1 = 0.0;
	sum2 = 0.0;
	invariant = 0.0;
	invariant_esd = 0.0;
	invariant_fit = 0.0;
	for (j = 1; j <= M; j++) {
		weight = pow( (double) dSAS[j], (double) -2 );
		sum1 += weight * (fit[j] - SAS[j]);
		sum2 += weight;
		if (j > 1) {
			Qsqr_dQ = Q[j] * Q[j] * (Q[j] - Q[j-1]) * 1.0e24;
			invariant     += Qsqr_dQ *  SAS[j];
			invariant_esd += Qsqr_dQ * dSAS[j];
			invariant_fit += Qsqr_dQ *  fit[j];
		}
	}
	suggested_background = background + sum1 / sum2;
	background_esd = 1.0 / sqrt(sum2);

#ifdef __MWERKS__
	fprintf (path, "\n");
	fprintf (path, "Summary statistics for size distribution\n");
	fprintf (path, "%15s\t %15s\t %15s\t %15s\n",
					"weighting", "number", "volume", "intensity" );
	fprintf (path, "%15s", "");
	if (vPower == 0) fprintf (path, "\t %15s", "chosen");
		else		 fprintf (path, "\t %15s", "derived");
	if (vPower == 1) fprintf (path, "\t %15s", "chosen");
		else		 fprintf (path, "\t %15s", "derived");
	if (vPower == 2) fprintf (path, "\t %15s", "chosen");
		else		 fprintf (path, "\t %15s", "derived");
	fprintf (path, "\n");
	fprintf (path, "\n");
	fprintf (path, "%15s\t %15g\t %15g\t %15g\n",
					"integral", nArea, vArea, iArea);
	fprintf (path, "%15s\t %15s\t %15s\t %15s\n",
					"", "1/cm^3", "fraction", "cm^3" );
	fprintf (path, "\n");
	fprintf (path, "%15s\t %15g\t %15g\t %15g\n",
					"mode, A", 2*nMode, 2*vMode, 2*iMode);
	fprintf (path, "%15s\t %15g\t %15g\t %15g\n", 
					"mean, A", 2*nMean, 2*vMean, 2*iMean);
	fprintf (path, "%15s\t %15g\t %15g\t %15g\n", 
					"std. dev., A", 2*nStdDev, 2*vStdDev, 2*iStdDev);
	if (AnalysisType == 1) {
	  fprintf (path, "\n");
	  fprintf (path, "%15s\t %15g\t %15g\t %15g\n", 
					"entropy", nEntropy, vEntropy, iEntropy);
	}

	fprintf (path, "\n");
	fprintf (path, "%40s\t %15g 1/cm\n", 
	    "background used for analysis", background);
	fprintf (path, "%40s\t %15g 1/cm\n", 
	    "suggested background", suggested_background);
	fprintf (path, "%40s\t %15g 1/cm\n", 
	    "estimated error of background", background_esd);
	fprintf (path, "%40s\n", 
	    "(based on analysis of weighted residuals)");

	fprintf (path, "\n");
	fprintf (path, "%40s\t %15g 1/cm^4\n", 
	    "SAS invariant calculated from fit", invariant_fit);
	fprintf (path, "%40s\t %15g 1/cm^4\n", 
	    "SAS invariant calculated from data", invariant);
	fprintf (path, "%40s\t %15g 1/cm^4\n", 
	    "esd of SAS invariant calculated from data", invariant_esd);
#else
	fprintf (path, "\n");
	fprintf (path, "Summary statistics for size distribution\n");
	fprintf (path, "%15s\t %15s\t %15s\t %15s\n",
					"weighting", "number", "volume", "intensity" );
	fprintf (path, "%15s", "");
	if (vPower == 0) fprintf (path, "\t %15s", "chosen");
		else		 fprintf (path, "\t %15s", "derived");
	if (vPower == 1) fprintf (path, "\t %15s", "chosen");
		else		 fprintf (path, "\t %15s", "derived");
	if (vPower == 2) fprintf (path, "\t %15s", "chosen");
		else		 fprintf (path, "\t %15s", "derived");
	fprintf (path, "\n");
	fprintf (path, "\n");
	fprintf (path, "%15s\t %15lg\t %15lg\t %15lg\n",
					"integral", nArea, vArea, iArea);
	fprintf (path, "%15s\t %15s\t %15s\t %15s\n",
					"", "1/cm^3", "fraction", "cm^3" );
	fprintf (path, "\n");
	fprintf (path, "%15s\t %15lg\t %15lg\t %15lg\n",
					"mode, A", 2*nMode, 2*vMode, 2*iMode);
	fprintf (path, "%15s\t %15lg\t %15lg\t %15lg\n", 
					"mean, A", 2*nMean, 2*vMean, 2*iMean);
	fprintf (path, "%15s\t %15lg\t %15lg\t %15lg\n", 
					"std. dev., A", 2*nStdDev, 2*vStdDev, 2*iStdDev);
	if (AnalysisType == 1) {
	  fprintf (path, "\n");
	  fprintf (path, "%15s\t %15lg\t %15lg\t %15lg\n", 
					"entropy", nEntropy, vEntropy, iEntropy);
	}

	fprintf (path, "\n");
	fprintf (path, "%40s\t %15lg 1/cm\n", 
	    "background used for analysis", background);
	fprintf (path, "%40s\t %15lg 1/cm\n", 
	    "suggested background", suggested_background);
	fprintf (path, "%40s\t %15lg 1/cm\n", 
	    "estimated error of background", background_esd);
	fprintf (path, "%40s\n", 
	    "(based on analysis of weighted residuals)");

	fprintf (path, "\n");
	fprintf (path, "%40s\t %15lg 1/cm^4\n", 
	    "SAS invariant calculated from fit", invariant_fit);
	fprintf (path, "%40s\t %15lg 1/cm^4\n", 
	    "SAS invariant calculated from data", invariant);
	fprintf (path, "%40s\t %15lg 1/cm^4\n", 
	    "esd of SAS invariant calculated from data", invariant_esd);
#endif 

       return (0);

}
