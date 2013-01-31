#include "utils.h"
#include <stdio.h>

const double SQRT_2PI = 2.5066282746310005024157652848110452530069867406099;
const double CORR_LEN = 0.1;
const double DAMP_RATE   = 0.5;
const double SQRT_C0 = 1.0;
const double VARIANCE = 0.5;
const double dt = 0.1;
const double CORR_TIME = 1.0;

void init()
{
	CORRL2 = CORR_LEN*CORR_LEN;
	SIGMA = sqrt(VARIANCE);
    TIME_RATE = dt / CORR_TIME;
    SQRT_TR = sqrt( TIME_RATE );
    SQRT_2TR = sqrt( 2*TIME_RATE );
}

double randGaussian()
{
	/*
	static double V1, V2, S;
	static int phase = 1;
	double X;

	if(phase) {
		do {
			double U1 = rand01();
			double U2 = rand01();

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
	*/
	
	double S = 0.0, U, V;
	
	while (S == 0.0 || S >= 1.0)
	{
		U = 2*rand1()-1;
		V = 2*rand1()-1;
		S = U*U + V*V;
	}
	return U * sqrt(-2.0 * log(S) / S);
}

void computeA(struct Complex** a, int nTerms, int mid)
{
#if TIMECORR_ON
    double amid_r = a[mid][mid].real;
#endif
	int i, j;
	for (i = 0; i < nTerms; i++)
	{
		for (j = mid; j < nTerms; j++)
		{
#if TIMECORR_ON
			a[i][j].real = a[i][j].real * (1 - TIME_RATE) 
                           + randGauss(0, SQRT_TR);
			a[i][j].imag = a[i][j].imag * (1 - TIME_RATE)
                           + randGauss(0, SQRT_TR);
/**			a[i][j].real = a[i][j].real * exp(- TIME_RATE) 
                           + randGauss(0, sqrt((1-exp(-2.0*TIME_RATE))/2.0));
			a[i][j].imag = a[i][j].imag * exp(- TIME_RATE)
                           + randGauss(0, sqrt((1-exp(-2.0*TIME_RATE))/2.0));**/
#else
            a[i][j].real = randGauss(1, INV_SQRT_2);
            a[i][j].imag = randGauss(0, INV_SQRT_2);
#endif
			a[2*mid-i][2*mid-j].real = a[i][j].real;
			a[2*mid-i][2*mid-j].imag = -a[i][j].imag;
			/*
			printf("i1 = %d, j1 = %d\ni2 = %d, j2 = %d\nIm a1 = %f\nIm a2 = %f\n", 
					i, j, 2*mid-i, 2*mid-j,
					a[i][j].imag, a[3*mid-i][2*mid-j].imag);
			*/
		}
	}
#if TIMECORR_ON
    a[mid][mid].real = amid_r * (1 - TIME_RATE)
                       + randGauss(0, SQRT_2TR);
/**    a[mid][mid].real = amid_r * exp(- TIME_RATE)
                       + randGauss(0, sqrt((1-exp(-2.0*TIME_RATE))/2.0));**/
#else
    a[mid][mid].real = randGauss(0, 1.0);                  
#endif
	a[mid][mid].imag = 0.0;
}

double scalarField(struct Vector2* pos, struct Complex** a, 
                    int nTerms, int mid, double sigma)
{
	double sum = 0.0;
	int i, j;
	 double sum2 = 0.0;
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
					exp(-k2*CORRL2/4.0);
			//sum += a[i][j].real;
			//printf("angle = %f    (ska vara 0)\n", angle);
			//printf("sum = %f    ==    %f = a.real  ??\n", sum*, a[i][j].real);
		}
	}
	//printf("sum = %f\n", sum);
	return SQRT_2PI * CORR_LEN * SQRT_C0 * sum;
}

void flowVel(struct Vector2* pos, struct Complex** a, 
        int nTerms, struct Vector2* u, int mid, double sigma)
{
    u->x = 0.0;
    u->y = 0.0;

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
            double ang = kx*pos->x + ky*pos->y;
            u->y += (-ky*a[i][j].real*sin(ang) - ky*a[i][j].imag*cos(ang) ) * exp(-k2*CORRL2/4.0);
            u->x +=  (-kx*a[i][j].real*sin(ang) - kx*a[i][j].imag*cos(ang) ) * exp(-k2*CORRL2/4.0);
        }
    }
    double coeff = SQRT_2PI * CORR_LEN * SQRT_C0;
    u->x *= coeff;
    u->y *= coeff;
    u->y = -u->y;
}

// A should be an array of length 4, [u_xx, u_xy; u_yx, u_yy]
void accMatrix(struct Vector2* pos, double A[], struct Complex** a, 
        int nTerms, int mid, double sigma)
{
	double sum = 0.0;
    A[0] = 0.0; A[1] = 0.0; A[2] = 0.0; A[3] = 0.0;
    double coeff = SQRT_2PI * CORR_LEN * SQRT_C0;
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
            double c = (a[i][j].real*cos(angle) - a[i][j].imag*sin(angle)) * 
					exp(-k2*CORRL2/4.0);
            A[0] += c * (-ky*kx);
            A[1] += c * (-ky*ky);
            A[2] += c * kx*kx;
            A[3] += c * kx*ky;
		}
	}
	A[0] *= coeff; A[1] *= coeff;
    A[2] *= coeff; A[3] *= coeff;
}

double factorial(int k)
{
    double fac = 1.0;
    int i;
    for (i = k; i > 1; i--)
    {
        fac *= i;
    }
    return fac;
}
