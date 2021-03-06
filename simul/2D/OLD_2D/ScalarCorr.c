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

	printf("\n*********************************************************\nTest correlation of velocity field.\n*********************************************************\n");
	
    FILE *file = fopen("velcorr", "w");
    
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

	printf("core dumped?");	        
	for (p2.x = -L/2.0; p2.x <= L/2.0; p2.x += dx)
	//for (p2.x = 0.0; p2.x <= L; p2.x += dx)
	{
		for (p2.y = -L/2.0; p2.y <= L/2.0; p2.y += dx)
		//for (p2.y = 0.0; p2.y <= L; p2.y += dx)
		{
			for (i = 0; i < n; i++)
			{
				computeA(a, nTerms, mid);
				struct Vector2 u1;
				struct Vector2 u2;
				flowVel(&p1, a, nTerms, &u1, mid, sigma);
				flowVel(&p2, a, nTerms, &u2, mid, sigma);
                sumUx += u1.y * u2.y / n;
                sumUy += u1.x * u2.x / n;
			}
			
            double Rx = p1.x - p2.x;
            double Ry = p1.y - p2.y;
			double expon = exp(-(Rx*Rx + Ry*Ry)/2.0/CORR_LENGTH2);
            double cffx = 1.0/CORR_LENGTH/CORR_LENGTH - 1.0/pow(CORR_LENGTH,4)*Ry*Ry;
            double cffy = 1.0/CORR_LENGTH/CORR_LENGTH - 1.0/pow(CORR_LENGTH,4)*Rx*Rx;

			    
			
            fprintf(file, "%.3f\t%.3f\t%.30f\t%.30f\t%.30f\t%.30f\n", 
                        p2.x, p2.y, sumUx, sumUy, expon*cffx, expon*cffy);
			printf("x = %f,  y = %f\n", p2.x, p2.y);
			printf("<psi*psi'> = %.10f\t", sumUx);
			printf("e... = %.10f\n\n", expon*cffx);
			sumUx = 0.0;
            sumUy = 0.0;
		}
		percent += dx / L;
		printf("%.2f %%\n", 100*percent);
	}
    fclose( file );
}

void testAccCorr(struct Complex** a, int nTerms, int mid, double sigma)
{

	printf("\n*********************************************************\nTest correlation of acc field.\n*********************************************************\n");
	
    FILE *file = fopen("acccorr", "w");
    
    double axx = 0.0, ayx = 0.0, axy = 0.0;
    double a12 = 0.0, a34 = 0.0, a14 = 0.0, a23 = 0.0;
	
	int n = 10000;
	
	struct Vector2 p1;
	struct Vector2 p2;
	
	p1.x = 0.0;
	p1.y = 0.0;
	
	int i;
	double x0 = 0.0, y0 = 0.0;
	
	double L = 1.0;
	double dx = 0.1;
	
	double percent = 0.0;

	for (p2.x = -L/2.0; p2.x <= L/2.0; p2.x += dx)
	//for (p2.x = 0.0; p2.x <= L; p2.x += dx)
	{
		for (p2.y = -L/2.0; p2.y <= L/2.0; p2.y += dx)
		//for (p2.y = 0.0; p2.y <= L; p2.y += dx)
		{
			for (i = 0; i < n; i++)
			{
				computeA(a, nTerms, mid);
                double A[] = { 0.0, 0.0, 0.0, 0.0 };
                double B[] = { 0.0, 0.0, 0.0, 0.0 };
                accMatrix(&p1, A, a, nTerms, mid, sigma);
				accMatrix(&p2, B, a, nTerms, mid, sigma);
                axx += A[0]*B[0] / n; ayx += A[2]*B[2] / n; axy += A[1]*B[1] / n;
                a12 += A[0]*B[1] / n; a34 += A[2]*B[3] / n; a14 += A[0]*B[3] / n; a23 += A[1]*B[2] / n;
			}
			
            double Rx = p1.x - p2.x;
            double Ry = p1.y - p2.y;
			double expon = exp(-(Rx*Rx + Ry*Ry)/2.0/CORR_LENGTH2);
            double cA12 = expon * Rx/CORR_LENGTH2*(Ry*Ry*Ry/pow(CORR_LENGTH,6) - 3.0*Ry/pow(CORR_LENGTH,4));
            double cA34 = expon * (-Ry/CORR_LENGTH2)*(-Rx*Rx*Rx/pow(CORR_LENGTH,6) + 3.0*Rx/pow(CORR_LENGTH,4));
            double cA14 = expon * (Ry*Ry/CORR_LENGTH2-1.0) * (1.0-Rx*Rx/CORR_LENGTH2) / CORR_LENGTH2 / CORR_LENGTH2;
            double cA23 = expon * (Rx*Rx/CORR_LENGTH2-1.0) * (1.0-Ry*Ry/CORR_LENGTH2) / CORR_LENGTH2 / CORR_LENGTH2;
            double cAxx = expon * (1 - Rx*Rx/CORR_LENGTH2)/CORR_LENGTH2;
            double cAxy = expon * (pow(Ry,4)/pow(CORR_LENGTH,4)-6.0*pow(Ry,2)/CORR_LENGTH2+720.0)/pow(CORR_LENGTH,4);
            double cAyx = expon * (pow(Rx,4)/pow(CORR_LENGTH,4)-6.0*pow(Rx,2)/CORR_LENGTH2+720.0)/pow(CORR_LENGTH,4);
			
            fprintf(file, "%.3f\t%.3f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\n", 
                        p2.x, p2.y, axx, cAxx, axy, cAxy, ayx, cAyx, a12, cA12, a34, cA34, a14, cA14, a23, cA23);
            axx = 0.0; axy = 0.0; ayx = 0.0;
            a12 = 0.0; a34 = 0.0; a14 = 0.0; a23 = 0.0;
		}
		percent += dx / L;
		printf("%.2f %%\n", 100*percent);
	}
    fclose( file );
}
