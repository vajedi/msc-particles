#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include "random.h"
#include "simath.h"

#define TIMECORR_ON 0

const double SQRT_2PI;
const double CORR_LEN;
double CORRL2;
const double DAMP_RATE;
const double SQRT_C0;
const double VARIANCE;
double SIGMA;

const double dt;
const double CORR_TIME;
double TIME_RATE;
double SQRT_TR;
double SQRT_2TR;

struct Complex;
struct Vector2;
struct Particle;

void init();
double randGaussian();
double factorial(int k);
void computeA(struct Complex** a, int nTerms, int mid);
double turbScalar(struct Vector2* pos, struct Complex** a, int nTerms, int mid, double sigma);
void turbFlow(struct Vector2* pos, struct Complex** a, 
	int nTerms, struct Vector2* u, int mid, double sigma);
void turbFlowGrad(struct Vector2* pos, double A[], struct Complex** a, 
        int nTerms, int mid, double sigma);

typedef struct
{
	double real;
	double imag;
} Complex;

typedef struct
{
	double x;
	double y;
} Vector2;

typedef struct
{
    Vector2 pos;
    Vector2 vel;
    Vector2 acc;
} Particle;

#endif
