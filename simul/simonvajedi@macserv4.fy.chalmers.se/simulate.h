#ifndef SIMULATE_H
#define SIMULATE_H

#include "utils.h"
#include "random.h"
#include <stdio.h>

/** Start with particle velocity 0. 
	Update the acceleration and so on...
**/

void initParticles(struct Particle* particles, int nParticles, double Lx, double Ly);
void updateParticles(struct Particle* particles, int nParticles, 
	struct Complex** a, struct Vector2*,double dt,int nTerms, int mid, double sigma, double dx);
void updatePos(struct Particle* p, double dt);
void updateVel(struct Particle* p, double dt);
void updateAcc(struct Particle* p, struct Complex** a, double dt, struct Vector2*,
			int nTerms, int mid, double sigma, double dx);

#endif
