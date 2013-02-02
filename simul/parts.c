#include <stdio.h>
#include "utils.h"
#include "testcorr.h"
#include "anorm.h"


const int nTerms = 20+1;//(int)(2 / CORR_LENGTH) + 1;//20 + 1;     nTerms = 2L/eta (see p. 30, Lic Thesis, K.G.)




int main(int argc, char* argv[])
{
	randomSeed();

    //printf("Correlation length = %.2f\n", CORR_LENGTH);
	
    //if ( nTerms % 2 == 0 )
    //    nTerms += 1;

	int mid = (nTerms - 1)/2;

	struct Complex** a = malloc( nTerms*sizeof(struct Complex*) );
	int i;
	for (i = 0; i < nTerms; i++)
	{
		a[i] = (struct Complex*) malloc( nTerms*sizeof(struct Complex) );
	}
		
	init();
	
	//simulateParticles(a, nTerms, mid, SIGMA, 0.1);
	//testANorm(a, nTerms, mid);
	//testScalarCorr(a, nTerms, mid, SIGMA);
    //testVelCorr(a, nTerms, mid, SIGMA);
    testAccCorr(a, nTerms, mid, SIGMA);
    //testTimeCorr(a, nTerms, mid, SIGMA);
	
    
    /*
	printf("mid = %d\n", mid);
	printf("nTerms = %d\n", nTerms);
	
	printf("%f + i%f\n%f + i%f\n%f + i%f\n", 
			a[mid+3][mid-1].real, a[mid+3][mid-1].imag, a[mid-3][mid+1].real, a[mid-3][mid+1].imag, 
			a[mid][mid].real, a[mid][mid].imag);
	*/
	
	for (i = 0; i < nTerms; i++)
	{
		free( a[i] );
	}
    free( a );
    
    
    
    return 0;
}
