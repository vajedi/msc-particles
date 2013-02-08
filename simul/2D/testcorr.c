#include "testcorr.h"

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
								(p1.y-p2.y)*(p1.y-p2.y))/2.0/CORRL2);
			
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
			double expon = exp(-(Rx*Rx + Ry*Ry)/2.0/CORRL2);
            double cffx = 1.0/CORR_LEN/CORR_LEN - 1.0/pow(CORR_LEN,4)*Ry*Ry;
            double cffy = 1.0/CORR_LEN/CORR_LEN - 1.0/pow(CORR_LEN,4)*Rx*Rx;

			    
			
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

    int use2dgrid = 1;
    FILE *file; 

    double axx = 0.0, ayx = 0.0, axy = 0.0;
    double a12 = 0.0, a34 = 0.0, a14 = 0.0, a23 = 0.0;

    int n = 10000;

    struct Vector2 p1;
    struct Vector2 p2;

    p1.x = 0.0;
    p1.y = 0.0;

    int i;
    double x0 = 0.0, y0 = 0.0;



    double percent = 0.0;
    if ( !use2dgrid )
    {
        double L = 0.4;
        double dx = 0.01;
        file = fopen("acccorr", "w");
        p2.x = 0.17; 
        double Rx = p1.x - p2.x;
        double Rx2 = Rx*Rx;
        for (p2.y = 0.011; p2.y <= L; p2.y += dx)
        {
            for (i = 0; i < n; i++)
            {
                computeA(a, nTerms, mid);
                double A[] = { 0.0, 0.0, 0.0, 0.0 };
                double B[] = { 0.0, 0.0, 0.0, 0.0 };
                accMatrix(&p1, A, a, nTerms, mid, sigma);
                accMatrix(&p2, B, a, nTerms, mid, sigma);
                axx += A[0]*B[0]; ayx += A[2]*B[2]; axy += A[1]*B[1];
                a12 += A[0]*B[1]; a34 += A[2]*B[3]; a14 += A[0]*B[3]; a23 += A[1]*B[2];
            }
            axx /= n; ayx /= n; axy /= n; a12 /= n; a34 /= n; a14 /= n; a23 /= n;
            double Ry = p1.y - p2.y;
            double Ry2 = Ry*Ry;
            double expon = exp(-(Rx2 + Ry2)/2.0/CORRL2);
            double cA12 = expon * Rx/CORRL2*(Ry2*Ry/pow(CORRL2,3) - 3.0*Ry/pow(CORR_LEN,4));
            double cA34 = expon * (-Ry/CORRL2)*(-Rx2*Rx/pow(CORR_LEN,6) + 3.0*Rx/pow(CORR_LEN,4));
            double cA14 = expon * (Ry2/CORRL2-1.0) * (1.0-Rx2/CORRL2) / CORRL2 / CORRL2;
            double cA23 = expon * (Rx2/CORRL2-1.0) * (1.0-Ry2/CORRL2) / CORRL2 / CORRL2;
            double cAxx = expon * (1.0 - Rx2/CORRL2) * (1.0 - Ry2/CORRL2) / CORRL2 / CORRL2;
            double cAxy = expon * (Ry2*Ry2/pow(CORR_LEN,4)-6.0*Ry2/CORRL2+3.0) / CORRL2 / CORRL2;
            double cAyx = expon * (Rx2*Rx2/pow(CORR_LEN,4)-6.0*Rx2/CORRL2+3.0) / CORRL2 / CORRL2;

            fprintf(file, "%.3f\t%.3f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\n", 
                    p2.x, p2.y, axx, cAxx, axy, cAxy, ayx, cAyx, a12, cA12, a34, cA34, a14, cA14, a23, cA23);
            axx = 0.0; axy = 0.0; ayx = 0.0;
            a12 = 0.0; a34 = 0.0; a14 = 0.0; a23 = 0.0;
            percent += dx / L;
            printf("%.2f %%\n", 100*percent);
        }
    }
    else
    {
    double L = 0.4;
    double dx = 0.04;
        file = fopen("acccorr2d","w");
        for (p2.x = -L/2.0; p2.x <= L/2.0; p2.x += dx)
        {
            double Rx = p1.x - p2.x;
            double Rx2 = Rx*Rx;
            for (p2.y = 0.011; p2.y <= L; p2.y += dx)
            {
                for (i = 0; i < n; i++)
                {
                    computeA(a, nTerms, mid);
                    double A[] = { 0.0, 0.0, 0.0, 0.0 };
                    double B[] = { 0.0, 0.0, 0.0, 0.0 };
                    accMatrix(&p1, A, a, nTerms, mid, sigma);
                    accMatrix(&p2, B, a, nTerms, mid, sigma);
                    axx += A[0]*B[0]; ayx += A[2]*B[2]; axy += A[1]*B[1];
                    a12 += A[0]*B[1]; a34 += A[2]*B[3]; a14 += A[0]*B[3]; a23 += A[1]*B[2];
                }
                axx /= n; ayx /= n; axy /= n; a12 /= n; a34 /= n; a14 /= n; a23 /= n;
                double Ry = p1.y - p2.y;
                double Ry2 = Ry*Ry;
                double expon = exp(-(Rx2 + Ry2)/2.0/CORRL2);
                double cA12 = expon * Rx/CORRL2*(Ry2*Ry/pow(CORRL2,3) - 3.0*Ry/pow(CORR_LEN,4));
                double cA34 = expon * (-Ry/CORRL2)*(-Rx2*Rx/pow(CORR_LEN,6) + 3.0*Rx/pow(CORR_LEN,4));
                double cA14 = expon * (Ry2/CORRL2-1.0) * (1.0-Rx2/CORRL2) / CORRL2 / CORRL2;
                double cA23 = expon * (Rx2/CORRL2-1.0) * (1.0-Ry2/CORRL2) / CORRL2 / CORRL2;
                double cAxx = expon * (1.0 - Rx2/CORRL2) * (1.0 - Ry2/CORRL2) / CORRL2 / CORRL2;
                double cAxy = expon * (Ry2*Ry2/pow(CORR_LEN,4)-6.0*Ry2/CORRL2+3.0) / CORRL2 / CORRL2;
                double cAyx = expon * (Rx2*Rx2/pow(CORR_LEN,4)-6.0*Rx2/CORRL2+3.0) / CORRL2 / CORRL2;

                fprintf(file, "%.3f\t%.3f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\t%.30f\n", 
                        p2.x, p2.y, axx, cAxx, axy, cAxy, ayx, cAyx, a12, cA12, a34, cA34, a14, cA14, a23, cA23);
                axx = 0.0; axy = 0.0; ayx = 0.0;
                a12 = 0.0; a34 = 0.0; a14 = 0.0; a23 = 0.0;
                
            }percent += dx / L;
                printf("%.2f %%\n", 100*percent);
        }
    }


    fclose( file );
}

/** Test time correlation using <a(t1)a(t2)> = exp(-|t2-t1|/tau) **/
void testTimeCorr(struct Complex** a, int nTerms, int mid, double sigma)
{
    printf("\n*********************************************************\nTest time correlation.\n*********************************************************\n");

    FILE *file = fopen("timecorr", "w");

    int numAvgs = 1000;
    int numTimeSteps = 50;

    computeA(a, nTerms, mid);

    int i, iAvg;
    struct Complex** a0 = malloc( nTerms*sizeof(struct Complex*) );
    double** sr = (double**) malloc( nTerms*sizeof(double*) );
    double** si = (double**) malloc( nTerms*sizeof(double*) );
    double* tvecr = (double*) malloc( numTimeSteps*sizeof(double) );
    double* tveci = (double*) malloc( numTimeSteps*sizeof(double) );
    double* et = (double*) malloc( numTimeSteps*sizeof(double) );
    for (i = 0; i < numTimeSteps; i++)
    {
        tvecr[i] = 0.0; tveci[i] = 0.0; et[i] = 0.0;
    }
    for (i = 0; i < nTerms; i++)
    {
        a0[i] = (struct Complex*) malloc( nTerms*sizeof(struct Complex) );
        sr[i]  = (double*) malloc( nTerms*sizeof(double*) );
        si[i]  = (double*) malloc( nTerms*sizeof(double*) );
        int j = 0;
        for (j = 0; j < nTerms; j++)
        {
            sr[i][j] = 0.0; si[i][j] = 0.0;
            a0[i][j].real = a[i][j].real;
            a0[i][j].imag = a[i][j].imag;
        }
    }

    int ix = 3, iy = 4;
    for (iAvg = 0; iAvg < numAvgs; iAvg++)
    {
        a0[ix][iy].real = a[ix][iy].real;
        a0[ix][iy].imag = a[ix][iy].imag;
        for (i = 0; i < numTimeSteps; i++)
        {
            //             for (ix = 0; ix < nTerms; ix++)
            //             {
            //                 for (iy = 0; iy < nTerms; iy++)
            //                 {
            tvecr[i] += (a[ix][iy].real*a0[ix][iy].real) / numAvgs;
            tveci[i] += (a[ix][iy].imag*a0[ix][iy].imag) / numAvgs;
            computeA(a, nTerms, mid);
            //                 }
            //             }
        }
        printf("%.1f %%\n", (double)(iAvg+1)*100.0/numAvgs);
    }


    for (i = 0; i < numTimeSteps; i++)
    {
        et[i] = exp(-dt*i/CORR_TIME) / 2.0;
        fprintf(file, "%.30f\t%.30f\t%.30f\n", tvecr[i], tveci[i], et[i]);
    }

    fclose( file );

    for (i = 0; i < nTerms; i++)
    {
        free( sr[i] );
        free( si[i] );
        free( a0[i] );
    }
    free( sr );
    free( si );
    free( tvecr );
    free( tveci );
    free( et );
    free( a0 );
}

/** Test time correlation using <a(t1)a(t2)> = exp(-|t2-t1|/tau) **/
void testTimeCorr2(struct Complex** a, int nTerms, int mid, double sigma)
{
    printf("\n*********************************************************\nTest time correlation.\n*********************************************************\n");

    FILE *file = fopen("timecorr2", "w");

    int numAvgs = 100000;
    int numTimeSteps = 50;

    struct Vector2 p1;
    struct Vector2 p2;

    double L = 0.5;
    double dx = L / (numTimeSteps - 1.0);
    p1.x = 0.0; p1.y = 0.0;
    p2.x = -L/2; p2.y = 0.2;
    p2.x = 0.05;
    int i, iAvg;


    double* svec = (double*) malloc( numTimeSteps*sizeof(double) );
    double* ux = (double*) malloc( numTimeSteps*sizeof(double) );
    double* uy = (double*) malloc( numTimeSteps*sizeof(double) );
    double CORRL4 = CORRL2*CORRL2;
    for (i = 0; i < numTimeSteps; i++)
    {
        svec[i] = 0; ux[i] = 0; uy[i] = 0;
    }

    for (iAvg = 0; iAvg < numAvgs; iAvg++)
    {    
        computeA(a, nTerms, mid);
        struct Vector2 u1;
        flowVel(&p1, a, nTerms, &u1, mid, sigma);
        double s1 = scalarField(&p1, a, nTerms, mid, sigma);

        for (i = 0; i < numTimeSteps; i++)
        {   
            double s2 = scalarField(&p2, a, nTerms, mid, sigma);
            svec[i] += s1*s2;
            struct Vector2 u2;
            flowVel(&p2, a, nTerms, &u2, mid, sigma);
            ux[i] += u1.x * u2.x;
            uy[i] += u1.y * u2.y;
            computeA(a, nTerms, mid);
        }
        //p2.x += dx;
        printf("%.1f %%\n", (double)(iAvg+1)*100.0/(numAvgs));
    }

    for (i = 0; i < numTimeSteps; i++)
    {
        svec[i] /= numAvgs; ux[i] /= numAvgs; uy[i] /= numAvgs;
        double Rx = p1.x - p2.x;
        double Ry = p1.y - p2.y;
        double expon = exp(-(Rx*Rx + Ry*Ry)/2.0/CORRL2-i*TIME_RATE);
        double cffx = 1.0/CORRL4 - 1.0/CORRL4*Ry*Ry;
        double cffy = 1.0/CORRL4 - 1.0/CORRL4*Rx*Rx;
        fprintf(file, "%.31f\t%.30f\t%.31f\t%.30f\t%.30f\t%.30f\n", svec[i], expon, ux[i], expon*cffx, uy[i], expon*cffy);
        //fprintf(file, "%.31f\t%.30f\t%.31f\t%.30f\t%.30f\t%.30f\n", svec[i], expon, sumUx, expon*cffx, sumUy, expon*cffy);
    }

    fclose( file );
    free( svec ); free( ux ); free( uy );
}
