#ifndef SCALARCORR_H
#define SCALARCORR_H

#include <stdio.h>
#include "utils.h"
#include "random.h"

void testScalarCorr(struct Complex** a, int nTerms, int mid, double sigma);

void testVelCorr(struct Complex** a, int nTerms, int mid, double sigma);
void testAccCorr(struct Complex** a, int nTerms, int mid, double sigma);
void testTimeCorr(struct Complex** a, int nTerms, int mid, double sigma);
#endif
