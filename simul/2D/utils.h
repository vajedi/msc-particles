#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include "random.h"
#include "simath.h"

const double SQRT_2PI;
const double CORR_LENGTH;
double CORR_LENGTH2;
const double DAMP_RATE;
const double SQRT_C0;
const double VARIANCE;
double SIGMA;

struct Complex;
struct Vector2;
struct Particle;

void init();
double randGaussian();
double factorial(int k);
void computeA(struct Complex** a, int nTerms, int mid);
double scalarField(struct Vector2* pos, struct Complex** a, int nTerms, int mid, double sigma);
void flowVel(struct Vector2* pos, struct Complex** a, 
	int nTerms, struct Vector2* u, int mid, double sigma);
void accMatrix(struct Vector2* pos, double A[], struct Complex** a, 
        int nTerms, int mid, double sigma);

struct Complex
{
	double real;
	double imag;
};

struct Vector2
{
	double x;
	double y;
};

struct Particle
{
		struct Vector2 pos;
		struct Vector2 vel;
		struct Vector2 acc;
};

#endif
