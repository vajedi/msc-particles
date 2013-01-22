#include "simulate.h"

//const double Dt = 0.1;
const int NBR_OF_PARTICLES = 100;
const double L = 1.0;

void simulateParticles(struct Complex** a, int nTerms, int mid, double sigma, double dx)
{
    printf("\n* * * * * * * * \n    Simulate the particles with Stokes' drag.\n* * * * * * * * *\n\n");
	struct Particle* particles = malloc(NBR_OF_PARTICLES*sizeof(struct Particle));
    struct Vector2* u = malloc( sizeof(struct Vector2) );
	initParticles(particles, NBR_OF_PARTICLES, L, L);
	double dt = 0.1;
	int i;
	for (i = 0; i < 100; i++)
	{
		updateParticles(particles, NBR_OF_PARTICLES, a, u, dt, nTerms, mid, sigma, dx);
	}
	
    printf("Done simulating particles. Freeing particles.\n\n");
	free(particles);
}

/** Distibute the particles randomly in the beginning **/
void initParticles(struct Particle* particles, int nParticles, double Lx, double Ly)
{
	int i;
	for (i = 0; i < nParticles; i++)
	{
		particles[i].pos.x = randd(-Lx/2, Lx/2);
		particles[i].pos.y = randd(-Ly/2, Ly/2);
		
		particles[i].vel.x = 0.0;
		particles[i].vel.y = 0.0;
		
		particles[i].acc.x = 0.0;
		particles[i].acc.y = 0.0;
	}
}

void updateParticles(struct Particle* particles, int nParticles, 
	struct Complex** a, struct Vector2* u, double dt,int nTerms, 
    int mid, double sigma, double dx)
{
	int i;
	for (i = 0; i < nParticles; i++)
	{
		particles[i].acc.x = 0;
		particles[i].acc.y = 0;

		updateAcc( &particles[i], a, dt, u, nTerms, mid, sigma, dx );
		updateVel( &particles[i], dt );
		updatePos( &particles[i], dt );
	}
}

void updatePos(struct Particle* p, double dt)
{
	p->pos.x += p->vel.x * dt;
	p->pos.y += p->vel.y * dt;
}

void updateVel(struct Particle* p, double dt)
{
	p->vel.x += p->acc.x * dt;
	p->vel.y += p->acc.y * dt;
}

void updateAcc(struct Particle* p, struct Complex** a, double dt, 
			struct Vector2* u, int nTerms, int mid, double sigma, double dx)
{
	flowVel(&(p->pos), a, nTerms, u, mid, sigma);
	
	p->acc.x += DAMP_RATE * (u->x - p->vel.x);
	p->acc.y += DAMP_RATE * (u->y - p->vel.y);
}
