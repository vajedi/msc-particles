#include "ScalarCorr.h"

void testScalarCorr(struct Complex** a, int nTerms, int mid, double sigma)
{

	printf("\n*********************************************************\nTest correlation of scalar field.\n*********************************************************\n");
	
    FILE *file = fopen("scalarcorr", "w");
    
	double sumScalar = 0.0;
	
	int n = 10000;
	
	struct Vector2 p1;
	struct Vector2 p2;
	
	p1.x = 0.0;
	p1.y = 0.0;
	
	int i;
	double x0 = 0.0, y0 = 0.0;
	
	double L = 1.0;
	double dx = 0.1;
	
	double percent = -dx/L;
	p2.x = 0.2;
	// for (p2.x = -L/2.0; p2.x <= L/2.0; p2.x += dx)
	//for (p2.x = 0.0; p2.x <= L; p2.x += dx)
	// {
		for (p2.y = -L/2.0; p2.y <= L/2.0; p2.y += dx)
		//for (p2.y = 0.0; p2.y <= L; p2.y += dx)
		{
			for (i = 0; i < n; i++)
			{
				computeA(a, nTerms, mid);
				double s1 = scalarField(&p1, a, nTerms, mid, sigma);
				double s2 = scalarField(&p2, a, nTerms, mid, sigma);
		
				sumScalar += s1*s2 / n;
			}
			double expon = exp(-((p1.x-p2.x)*(p1.x-p2.x)+
								(p1.y-p2.y)*(p1.y-p2.y))/2.0/CORR_LENGTH2);
			
            fprintf(file, "%.3f\t%.3f\t%.15f\t%.15f\n", 
                        p2.x, p2.y, sumScalar, expon);
			printf("x = %f,  y = %f\n", p2.x, p2.y);
			printf("<psi*psi'> = %.10f\t", sumScalar);
			printf("e... = %.10f\n\n", expon);
			sumScalar = 0.0;
		}
		percent += dx / L;
		printf("%.2f %%\n", 100*percent);
	// }
    fclose( file );
}


// Se ekv. (4.3)
void testVelCorr(struct Complex** a, int nTerms, int mid, double sigma)
{

	printf("\n*********************************************************\nTest correlation of scalar field.\n*********************************************************\n");
	
    FILE *file = fopen("scalarcorr", "w");
    
    double sumUx = 0.0; // Test u_x (d_y phi)
    double sumUy = 0.0;
	
	int n = 10000;
	
	struct Vector2 p1;
	struct Vector2 p2;
	
	p1.x = 0.0;
	p1.y = 0.0;
	
	int i;
	double x0 = 0.0, y0 = 0.0;
	
	double L = 1.0;
	double dx = 0.1;
	
	double percent = -dx/L;
	p2.x = 0.2;
	for (p2.x = -L/2.0; p2.x <= L/2.0; p2.x += dx)
	//for (p2.x = 0.0; p2.x <= L; p2.x += dx)
	{
		for (p2.y = -L/2.0; p2.y <= L/2.0; p2.y += dx)
		//for (p2.y = 0.0; p2.y <= L; p2.y += dx)
		{
			for (i = 0; i < n; i++)
			{
				computeA(a, nTerms, mid);
				struct Vector2* u1;
				struct Vector2* u2;
				flowVel(&p1, a, nTerms, u1, mid, sigma);
				flowVel(&p2, a, nTerms, u2, mid, sigma);
		        
                sumUx += u1.x * u2.x / n;
                sumUy += u1.y * u2.y / n;
			}
			
			double expon = exp(-((p1.x-p2.x)*(p1.x-p2.x)+
								(p1.y-p2.y)*(p1.y-p2.y))/2.0/CORR_LENGTH2);

			for (i = 0; i < 2; i++) // 2 because of 2 dimensions
			{
			    //double psum = 0.0;
			    //int pind;
			    //for (pind = 0; pind < 
			    
			
            fprintf(file, "%.3f\t%.3f\t%.15f\t%.15f\n", 
                        p2.x, p2.y, sumUx, sumUy, expon);
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
