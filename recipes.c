/* recipes.c

    Routines from Numerical Recipes
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "recipes.h"

double  *vector(int nl, int nh);
void    free_vector (double *v, int nl, int nh);
int     *ivector(int nl, int nh);
void    free_ivector (int *v, int nl, int nh);
double  **matrix(int nrl, int nrh, int ncl, int nch);
void    free_matrix (double **m, int nrl, int nrh, int ncl, int nch);
void    hunt(double xx[], int n, double x, int *jlo);
void    nrerror (char error_text[]);


double *vector(int nl, int nh)  /* allocate a vector double[nl..nh] */
{
    double *v;
    v = (double *) malloc ((unsigned long) (nh-nl+1)*sizeof(double));
    if (!v) nrerror ("memory allocation failure in vector()");
    return v-nl;
}

void free_vector (double *v, int nl, int nh)
{
    free ((char*) (v+nl));
}

int *ivector(int nl, int nh)    /* allocate a vector int[nl..nh] */
{
    int *v;
    v = (int *) malloc ((unsigned long) (nh-nl+1)*sizeof(int));
    if (!v) nrerror ("memory allocation failure in ivector()");
    return v-nl;
}

void free_ivector (int *v, int nl, int nh)
{
    free ((char*) (v+nl));
}

double **matrix(int nrl, int nrh, int ncl, int nch)
{
    int i, qty;
    double **m;
    m = (double **) malloc ((unsigned long) (nrh-nrl+1)*sizeof(double*));
    if (!m) nrerror ("memory allocation failure #1 in matrix()");
    m -= nrl;
    qty = nch-ncl+1;
    for (i=nrl;i<=nrh;i++) {
        m[i] = (double *) malloc ((unsigned long) (qty)*sizeof(double));
        if (!m[i]) nrerror ("memory allocation failure #2 in matrix()");
        m[i] -= ncl;
    }
    return m;
}

void free_matrix (double **m, int nrl, int nrh, int ncl, int nch)
{
    int i;
    for (i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
    free((char*) (m+nrl));
}

/************
 *** hunt ***       Numerical Recipes routine
 ************
    Find xx[j] < x <= xx[j+1] where 1 <= jlo <= n.  
      If x <= xx[1], jlo = 1.
      If x > xx[n], jlo = n + 1.
      Returns the value of jlo.
    x: value to search
    xx: array of table values, numbered 1 to n
    jlo: returned index (also used as an input)
*/
void hunt(double xx[], int n, double x, int *jlo)
{
    int jm,jhi,inc,ascnd;
    ascnd=(xx[n] > xx[1]);
    if (*jlo <= 0 || *jlo > n) {
        *jlo = 0;
        jhi = n+1;
    } else {
        inc = 1;
        if (x >= xx[*jlo] == ascnd) {
            if (*jlo == n) return;
            jhi = (*jlo)+1;
            while (x >= xx[jhi] == ascnd) {
                *jlo = jhi;
                inc += inc;
                jhi = (*jlo)+inc;
                if (jhi > n) {
                    jhi = n+1;
                    break;
                }
            }
        } else {
            if (*jlo ==1) {
                *jlo = 0;
                return;
            }
            jhi = (*jlo);
            *jlo -= 1;
            while (x < xx[*jlo] == ascnd) {
                jhi = (*jlo);
                inc += inc;
                *jlo = jhi - inc;
                if (*jlo < 1) {
                    *jlo = 0;
                    break;
                }
            }
        }
    }
    while (jhi-(*jlo) != 1) {
        jm = (jhi+(*jlo)) >> 1;
        if (x > xx[jm] == ascnd)
            *jlo = jm;
        else
            jhi = jm;
    }
}

void nrerror (char error_text[])        /* NR standard error handler */
{
    fprintf (stderr, "Numerical Recipes run-time error...\n");
    fprintf (stderr, "%s\n", error_text);
    fprintf (stderr, "... now exiting to system ...\n");
    /*  WaitMouseClick ();  */
    exit (1);   /* return a non-zero value */
}

/**********************************************************************/
/*  svdcmp.c */

  /* 2004-12-23,PRJ:
    guard against underflow of g leading to INF or NAN
   */
#define G_UNDERFLOW  (1.0e-280)

void svdcmp(double **a, int m, int n, double w[], double **v)
{
	double pythag(double a, double b);
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=vector(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
        if (fabs(g)<G_UNDERFLOW) g = 0;
        for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];  if (fabs(g)<G_UNDERFLOW) g = 0;
		l=i;
	}
	for (i=MIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];  if (fabs(g)<G_UNDERFLOW) g = 0;
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;  if (fabs(g)<G_UNDERFLOW) g = 0;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];  if (fabs(g)<G_UNDERFLOW) g = 0;
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];  if (fabs(g)<G_UNDERFLOW) g = 0;
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);  if (fabs(g)<G_UNDERFLOW) g = 0;
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];  if (fabs(g)<G_UNDERFLOW) g = 0;
				y=w[i];
				h=s*g;
				g=c*g;  if (fabs(g)<G_UNDERFLOW) g = 0;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1,1,n);
}


/**********************************************************************/
/*  svbksb.c */

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
	int jj,j,i;
	double s,*tmp;

	tmp=vector(1,n);
	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
if (isnan(s)) {
  printf("<%s,%d> x[%d]=%lg\n", __FILE__, __LINE__, j, s);
  for (jj=1;jj<=n;jj++) 
    printf (
      "<%s,%d> (v[%d][%d]=%lg)*(tmp[%d]=%lg) = %lg\n", 
      __FILE__, __LINE__, j, jj, v[j][jj], jj, tmp[jj], v[j][jj]*tmp[jj]
    );
    exit(0);
}
	}
	free_vector(tmp,1,n);
}


/**********************************************************************/
/*  pythag.c */

double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


/**********************************************************************/
/*  gasdev.c */

double gasdev(long *idum)
{
	double ran1(long *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

/**********************************************************************/
/*  ran1.c */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/**********************************************************************/
