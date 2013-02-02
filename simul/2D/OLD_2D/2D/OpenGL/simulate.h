#ifndef SIMULATE_H
#define SIMULATE_H

#include "utils.h"
#include "random.h"
#include "simath.h"
#include <GL/glut.h>
/** Start with particle velocity 0. 
	Update the acceleration and so on...
**/

struct Particle* particles;
struct Vector2* u;

const int NBR_OF_PARTICLES;
const double L;

void initParticles(int, int, double);
void update();
void updateParticles(struct Particle* particles, int nParticles, 
                        struct Vector2* u);
void updatePos(struct Particle* p, double dt);
void updateVel(struct Particle* p, double dt);
void updateAcc(struct Particle* p, double dt, struct Vector2*,
			int nTerms, int mid, double sigma, double dx);
void draw();

#endif
