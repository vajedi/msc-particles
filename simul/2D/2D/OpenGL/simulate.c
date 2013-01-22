#include "simulate.h"
#include <stdio.h>

int Mid = 0;
int NTerms = 0;
double sigma;

void update()
{
    computeA(NTerms, Mid);
    updateParticles(particles, NBR_OF_PARTICLES, u);
    
    glutPostRedisplay();
    glutTimerFunc(1, update, 0);
}

float maxSpeed = 0.0;
float scalex = 1.5f;
float scaley = 1.5f;

void draw()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glBegin(GL_POINTS);
    int i;
    for (i = 0; i < NBR_OF_PARTICLES; i += 10)
    {
        float speed = (float)lengthSquared( &particles[i].vel );
        if (speed > maxSpeed)
            maxSpeed = speed;
    }
    
    for (i = 0; i < NBR_OF_PARTICLES; i++)
    {
        float speed = (float)lengthSquared( &particles[i].vel );
        float q = speed/maxSpeed;
        float q2 = q*q;
        //glColor3f(max(0,-4.0f*q2+8.0f*q-3.0f), max(0,4.0f*q-4.0f*q2), //8.0f*q-8.0f*q2-1.5f    //4.0f*q-4.0f*q2
        //          max(0,1.0f-8.0f*q2));
        glColor3f(min(1,q), max(0,1.0f-q), 0);
        glVertex3f(scalex*particles[i].pos.x, scaley*particles[i].pos.y, 0);
    }
    glEnd();
    
    glutSwapBuffers();
}

/** Distibute the particles randomly in the beginning **/
void initParticles(int _nTerms, int _mid, double s)
{
    NTerms = _nTerms;
    Mid = _mid;
    sigma = s;
    particles = malloc(NBR_OF_PARTICLES*sizeof(struct Particle));
    u = malloc( sizeof(struct Vector2) );
	int i;
	for (i = 0; i < NBR_OF_PARTICLES; i++)
	{
		particles[i].pos.x = randd(-L/2.0, L/2.0);
		particles[i].pos.y = randd(-L/2.0, L/2.0);
		
		particles[i].vel.x = 0.0;
		particles[i].vel.y = 0.0;
		
		particles[i].acc.x = 0.0;
		particles[i].acc.y = 0.0;
	}
}

void updateParticles(struct Particle* particles, int nParticles, 
	struct Vector2* u)
{
	int i;
	for (i = 0; i < nParticles; i++)
	{
		updatePos( &particles[i], dt );
		updateVel( &particles[i], dt );
		particles[i].acc.x = 0;
		particles[i].acc.y = 0;
		updateAcc( &particles[i], dt, u, NTerms, Mid, sigma, dx );
	}
}

void updatePos(struct Particle* p, double dt)
{
	p->pos.x += p->vel.x * dt;
	p->pos.y += p->vel.y * dt;
    
    p->pos.x = wrap( p->pos.x, -L/2.0, L/2.0);
    p->pos.y = wrap( p->pos.y, -L/2.0, L/2.0);
}

void updateVel(struct Particle* p, double dt)
{
	p->vel.x += p->acc.x * dt;
	p->vel.y += p->acc.y * dt;
}

void updateAcc(struct Particle* p, double dt, 
			struct Vector2* u, int nTerms, int mid, double sigma, double dx)
{
    /** Implement the Maxey-Riley equation here **/
	flowVel(&(p->pos), NTerms, u, Mid, sigma, dx);
	
	p->acc.x += DAMP_RATE * (u->x - p->vel.x);
	p->acc.y += DAMP_RATE * (u->y - p->vel.y);
}
