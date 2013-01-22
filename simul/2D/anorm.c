#include "anorm.h"

void testANorm(struct Complex** a, int nTerms, int mid)
{
	printf("Test if a is normalized (< a_k a_k* > = d_kk' \?).\n");
	int n = 10000, i, k;
	double sum = 0.0;
	for (i = 0; i < n; i++)
	{
		computeA(a, nTerms, mid);
		for (k = 0; k < mid; k++)
		{
			/*sum += (a[k][2*mid-k].real*a[k][2*mid-k].real + 
					a[k][2*mid-k].imag*a[k][2*mid-k].imag) / n / mid;*/
			sum += (a[k][k].real*a[1][2].real + 
					a[k][k].imag*a[1][2].imag) / n / mid;
		}
	}
	printf("%f\n", sum);
}