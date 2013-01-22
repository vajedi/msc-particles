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

const int NBR_OF_PARTICLES;
const double L;
double dt;
double dx;

struct Complex;
struct Vector2;
struct Particle;

struct Complex** a;
struct Complex** b;

void init(int, int);
void computeA(int nTerms, int mid);
double computeB(int nTerms, int mid);
double computeCoeffs(struct Complex** c, int nTerms, int mid, double stddev);
double scalarField(struct Vector2* pos, int nTerms, int mid, double sigma);
double derivScalarField(struct Vector2* pos, 
                    int nTerms, int mid, double sigma, int ind);
void flowVel(struct Vector2* pos, 
	int nTerms, struct Vector2* u, int mid, double sigma, double dx);
double lengthSquared(struct Vector2* v);
void freeA(int);

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
