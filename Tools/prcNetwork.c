
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void trapz_twoCell(double* Q, double* K, const unsigned N, const double dx, double* result)
{
	unsigned i, j, ij;
	double sum;
	double F[N*N];

	// set up the integrand
	for(i=0; i<N; i++)		// phi_1
	{
		for(j=0; j<N; j++)	// phi_2
		{
			F[i*N+j] = Q[i]*K[N*i+j]-Q[j]*K[N*j+i];	// F(phi_1, phi_2)
		}
	}

	// trapz integration
	for(i=0; i<N; i++)	// dphi index
	{
		for(j=0; j<N; j++)	// phi index (integration variable)
		{
			ij = fmod(j-i, N);	// phi_2 = dphi+phi
			sum = ( (double)j*sum + F[N*j+ij] )/(double)(j+1);
		}

		result[i] = dx*sum;
	}
}


void memoryError()
{
	printf("prcNetwork: Out of memory.");
	exit(-1);
}


void trapz_threeCell(double* Q, double* K, const unsigned N, const double dx, double* result)
{
	unsigned i, j, k, i12, i13, k12, k13, N2=N*N;
	double F1, sum12, sum13;
	double *F12, *F13;

	if((F12 = (double*)malloc(N*N2*sizeof(double))) == NULL) memoryError();
	if((F13 = (double*)malloc(N*N2*sizeof(double))) == NULL) memoryError();

	// set up the integrands
	for(i=0; i<N; i++)			// phi_1
	{
		for(j=0; j<N; j++)		// phi_2
		{
			for(k=0; k<N; k++)	// phi_3
			{
				F1 = Q[i]*(K[N*i+j]+K[N*i+k]);
				F12[i*N2+j*N+k] = F1-Q[j]*(K[N*j+i]+K[N*j+k]);	// F_12(phi_i, phi_j, phi_k)
				F13[i*N2+j*N+k] = F1-Q[k]*(K[N*k+i]+K[N*k+j]);	// F_13(phi_i, phi_j, phi_k)
			}
		}
	}

	// trapz integration
	for(i12=0; i12<N; i12++)		// dphi_12 index
	{
		for(i13=0; i13<N; i13++)	// dphi_13 index
		{
			for(k=0; k<N; k++)	// phi index (integration variable)
			{
				k12 = fmod(k-i12, N);	// phi - dphi12
				k13 = fmod(k-i13, N);	// phi - dphi13
				sum12 = ((double)k*sum12 + F12[k*N2+k12*N+k13]) / (double)(k+1);
				sum13 = ((double)k*sum13 + F13[k*N2+k12*N+k13]) / (double)(k+1);
			}
	
			result[i12*N+i13]    = dx*sum12;
			result[i12*N+i13+N2] = dx*sum13;
		}
	}
	free(F12);
	free(F13);
}
