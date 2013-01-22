#include "ScalarCorr.h"

void testScalarCorr(int nTerms, int mid, double sigma)
{

	printf("\n*********************************************************\nTest correlation of scalar field.\n*********************************************************\n");
	
    FILE *file = fopen("scalarcorr", "w");
    
	double sumScalar = 0.0;
	
	int n = 1000;
	
	struct Vector2 p1;
	struct Vector2 p2;
	
	p1.x = 0.0;
	p1.y = 0.0;
	
	int i;
	double x0 = 0.0, y0 = 0.0;
	
	double L = 1.0;
	double dx = 0.1;
	
	double percent = -dx/L;
	
	//for (p2.x = -L/2.0; p2.x <= L/2.0; p2.x += dx)
	for (p2.x = 0.0; p2.x <= L; p2.x += dx)
	{
		//for (p2.y = -L/2.0; p2.y <= L/2.0; p2.y += dx)
		for (p2.y = 0.0; p2.y <= L; p2.y += dx)
		{
			for (i = 0; i < n; i++)
			{
				computeA(nTerms, mid);
				double s1 = scalarField(&p1, nTerms, mid, sigma);
				double s2 = scalarField(&p2, nTerms, mid, sigma);
		
				sumScalar += s1*s2 / n;
			}
			double expon = VARIANCE*exp(-((p1.x-p2.x)*(p1.x-p2.x)+
								(p1.y-p2.y)*(p1.y-p2.y))/2.0/CORR_LENGTH2);
			
            fprintf(file, "%.3f\t%.3f\t%.15f\t%.15f\n", 
                        p2.x, p2.y, VARIANCE*sumScalar, expon);
			printf("x = %f,  y = %f\n", p2.x, p2.y);
			printf("<psi*psi'> = %.10f\t", sumScalar);
			printf("e... = %.10f\n\n", expon);
			sumScalar = 0.0;
		}
		percent += dx / L;
		printf("%.2f %%\n", 100*percent);
	}
    fclose( file );
}
