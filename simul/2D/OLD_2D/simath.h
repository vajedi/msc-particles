#ifndef SIMATH_H
#define SIMATH_H

/**
        Constants
\** **/
const double SQRT_2;
const double INV_SQRT_2;
const double PI;
const double TWO_PI;
const double PI_OVER_2;
const double PI_OVER_4;
const double PI_OVER_6;

/** Wrap a value periodically within a given range.
 * 
 * value - The variable to wrap.
 * min - The minimum value to return (inclusive).
 * max - The upper bound of the range (exclusive).
 **/
double wrap(double value, double min, double max);

/** Clamp a value within a given range.
 * 
 * value - The variable to clamp.
 * min - The minimum value of the range (inclusive).
 * max - The maximum value of the range (exclusive).
 **/
double clamp(double value, double min, double max);

/** Compute the linear interpolation between two values.
 * 
 * min - The minimum return value (inclusive).
 * max - The maximum return value (inclusive).
 * weight - A value in range [0,1] which is used 
 *          to interpolate between min and max.
 **/
double lerp(double min, double max, double weight);

#endif
