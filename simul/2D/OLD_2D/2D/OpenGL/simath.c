
const double SQRT_2 = 1.4142135623730950488016887242096980785696718753769480;
const double INV_SQRT_2 = 0.7071067811865475244008443621048490392848359376884740;
const double PI = 3.1415926535897932384626433832795028841971693993751058;
const double TWO_PI = 6.2831853071795864769252867665590057683943387987502116;
const double PI_OVER_2 = 1.5707963267948966192313216916397514420985846996875529;
const double PI_OVER_4 = 0.7853981633974483096156608458198757210492923498437764;
const double PI_OVER_6 = 0.5235987755982988730771072305465838140328615665625176;

double wrap(double value, double min, double max)
{
    double diff = max - min;
    if (diff < 0)
        return value;
    while (value < min)
        value += diff;
    while (value > max)
        value -= diff;
    return value;
}

double clamp(double value, double min, double max)
{
    if (min > max)
        return value;
    if (value < min)
        return min;
    if (value > max)
        return max;
    return value;
}

double lerp(double min, double max, double weight)
{
    return min + (max - min) * weight;
}

double max(double a, double b)
{
    if (a > b)
        return a;
    return b;
}

double min(double a, double b)
{
    if (a < b)
        return a;
    return b;
}
