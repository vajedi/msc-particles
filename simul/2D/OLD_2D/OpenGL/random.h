#ifndef RANDOM_H
#define RANDOM_H

#include <stdlib.h>
#include <math.h>
#include <time.h>

/** Use random seed. **/
void randomSeed();

/** Returns a random number wihin a given range.
 * 
 * min - The minimum value (inclusive).
 * max - The maximum value (exclusive).
**/
double randd(double min, double max);


/** Returns a random integer between.
 * 
 * min - The minimum value (inclusive).
 * max - The maximum value (exclusive).
**/
int randi(int min, int max);

/** Returns a random number between 0 and 1.
 * 
 * 0 is inclusive, and 1 exlusive.
**/
double rand1();

/** Return a Gaussian normal distributed random number.
**/
double randNormal();

/** Return a Gaussian distributed random number.
 *  
 * mean - The mean value of the Gaussian distribution.
 * dev  - The standard deviation of the Gaussian distribution.
**/
double randGauss(double mean, double dev);

#endif
