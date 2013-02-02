#include "random.h"

void randomSeed()
{
	srand((unsigned)time(NULL));
}

double randd(double min, double max)
{
	return (rand()/((double)RAND_MAX + 1)) * (max - min) + min;
}

int randi(int min, int max)
{
    return (int)(rand()/((double)RAND_MAX + 1)) 
                * (max - min) + min;
}

double rand1()
{
	return randd(0.0, 1.0);
}

double randNormal()
{
    double U, V, S = 0.0;
    while (S == 0 || S >= 1)
    {
        U = randd(-1, 1);
        V = randd(-1, 1);
        S = U*U + V*V;
    }
    return U * sqrt(-2.0 * log(S) / S);
}

double randGauss(double mean, double dev)
{
    double G = randNormal();
    return G*dev + mean;
}