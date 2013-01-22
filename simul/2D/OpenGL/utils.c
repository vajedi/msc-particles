#include "utils.h"
#include <stdio.h>

const double SQRT_2PI = 2.5066282746310005024157652848110452530069867406099;
const double CORR_LENGTH = 0.060;
const double CORR_TIME = 0.25;
const double DAMP_RATE   = 0.050;
const double SQRT_C0 = 1.0;
const double VARIANCE = 0.119;
double SQRT_dT_OVER_CORRTIME;
double TIME_DECAY;

const int NBR_OF_PARTICLES = 500;
const double L = 1.0;
double dt = 0.01;
double dx = 0.01;


struct Complex** a;
struct Complex** b;

void init(int nTerms, int mid)
{
    a = malloc( nTerms*sizeof(struct Complex*) );
    b = malloc( nTerms*sizeof(struct Complex*) );
	int i;
	for (i = 0; i < nTerms; i++)
	{
		a[i] = malloc( nTerms*sizeof(struct Complex) );
        b[i] = malloc( nTerms*sizeof(struct Complex) );
	}
    computeCoeffs(a, nTerms, mid, INV_SQRT_2);
	CORR_LENGTH2 = CORR_LENGTH*CORR_LENGTH;
	SIGMA = sqrt(VARIANCE);
    SQRT_dT_OVER_CORRTIME = sqrt(dt / CORR_TIME);
    TIME_DECAY = 1.0 - dt / CORR_TIME;
}

void computeA(int nTerms, int mid)
{
    computeB(nTerms, mid);
	int i, j;
	for (i = 0; i < nTerms; i++)
	{
		for (j = mid; j < nTerms; j++)
		{
			a[i][j].real = a[i][j].real*TIME_DECAY + b[i][j].real;
			a[i][j].imag = a[i][j].imag*TIME_DECAY + b[i][j].imag;
			
			a[2*mid-i][2*mid-j].real = a[i][j].real;
			a[2*mid-i][2*mid-j].imag = -a[i][j].imag;
		}
	}
	a[mid][mid].imag = 0.0;
}

double computeB(int nTerms, int mid)
{
    computeCoeffs(b, nTerms, mid, SQRT_dT_OVER_CORRTIME);
}

double computeCoeffs(struct Complex** c, int nTerms, int mid, double stddev)
{
	int i, j;
	for (i = 0; i < nTerms; i++)
	{
		for (j = mid; j < nTerms; j++)
		{
			c[i][j].real = randGauss(0, stddev);
			c[i][j].imag = randGauss(0, stddev);
			
			c[2*mid-i][2*mid-j].real = c[i][j].real;
			c[2*mid-i][2*mid-j].imag = -c[i][j].imag;
		}
	}
    c[mid][mid].real = randGauss(0, 2*stddev*stddev);
	c[mid][mid].imag = 0.0;
}

double scalarField(struct Vector2* pos, 
                    int nTerms, int mid, double sigma)
{
	double sum = 0.0;
	int i, j;
	for (i = 0; i < nTerms; i++)
	{
		double kx = i - mid;
		kx *= TWO_PI;
		for (j = 0; j < nTerms; j++)
		{
			double ky = j - mid;
			ky *= TWO_PI;
			double k2 = kx*kx + ky*ky;
			double angle = kx*pos->x + ky*pos->y;
			sum += (a[i][j].real*cos(angle) - a[i][j].imag*sin(angle)) * 
					exp(-k2*CORR_LENGTH2/4.0);
		}
	}
	return SQRT_2PI * CORR_LENGTH * SQRT_C0 * sum;
}

/** ind declares with respect to what variable the scalar field  
 *  should be differentiated. 0 if along x axis, 1 for y axis.
**/
double derivScalarField(struct Vector2* pos, 
                    int nTerms, int mid, double sigma, int ind)
{
	double sum = 0.0;
	int i, j;
	for (i = 0; i < nTerms; i++)
	{
        double k[] = { i - mid, 0 };
		k[0] *= TWO_PI;
		for (j = 0; j < nTerms; j++)
		{
			k[1] = j - mid;
			k[1] *= TWO_PI;
			double k2 = k[0]*k[0] + k[1]*k[1];
			double angle = k[0]*pos->x + k[1]*pos->y;
			sum += k[ind]*(-a[i][j].imag*cos(angle) - a[i][j].real*sin(angle)) * 
					exp(-k2*CORR_LENGTH2/4.0);
		}
	}
	return SQRT_2PI * CORR_LENGTH * SQRT_C0 * sum;
}

void flowVel(struct Vector2* pos, 
        int nTerms, struct Vector2* u, int mid, double sigma, double dx)
{
	u->x = derivScalarField(pos, nTerms, mid, sigma, 0);
	u->y = derivScalarField(pos, nTerms, mid, sigma, 1);
}

double lengthSquared(struct Vector2* v)
{
    return v->x * v->x + v->y * v->y;
}

void freeA(int nTerms)
{
    int i;
	for (i = 0; i < nTerms; i++)
	{
		free( a[i] );
        free( b[i] );
	}
    free( a );
    free( b );
}
