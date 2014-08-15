
/*
 */


#include <math.h>
#include <stdio.h>


#define N_EQ1		1
#define NP1		2
#define N_EQ3		3*N_EQ1
#define NP3		3*NP1
#define N_EQ4		4*N_EQ1
#define NP4		4*NP1
#define NC		6


#define THRESHOLD_SLOPE		10. //

extern "C" {

double boltzmann_sin(const double phase, const double k)
{
	return 1./(1.+exp(k*sin(phase)));
}

double boltzmann_cos(const double phase, const double k)
{
	return 1./(1.+exp(k*cos(phase)));
}


// p[0] = omega
// p[1] = beta
void derivs_one(double* y, double* dxdt, const double* p)
{
	dxdt[0] = p[0]-p[1]*cos(y[0])-cos(2.*y[0]);
}


void integrate_one_rk4(double* y, const double* params, double* output, const double dt, const unsigned N, const unsigned stride)
{
	unsigned i, j, k;
	double dt2, dt6;
	double y1[N_EQ1], y2[N_EQ1], k1[N_EQ1], k2[N_EQ1], k3[N_EQ1], k4[N_EQ1];
	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<N_EQ1; j++) output[j] = y[j];
	for(i=1; i<N; i++) {
		for(j=0; j<stride; j++) {
			derivs_one(y, k1, params);
			for(k=0; k<N_EQ1; k++) y1[k] = y[k]+k1[k]*dt2; 			
			derivs_one(y1, k2, params);
			for(k=0; k<N_EQ1; k++) y2[k] = y[k]+k2[k]*dt2; 			
			derivs_one(y2, k3, params);
			for(k=0; k<N_EQ1; k++) y2[k] = y[k]+k3[k]*dt; 			
			derivs_one(y2, k4, params);
			for(k=0; k<N_EQ1; k++)
			{
				y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
				y[k] = fmod(y[k], 2.*M_PI);
			}
		}
		for(j=0; j<N_EQ1; j++) output[N_EQ1*i+j] = y[j];
	}
}


// kij[0] :  1 -o  0
// kij[1] :  2 -o  0
// kij[2] :  0 -o  1
// kij[3] :  2 -o  1
// kij[4] :  0 -o  2
// kij[5] :  1 -o  2
void derivs_three(double* y, double* dxdt, const double* p, const double* kij)
{
	double bm_cos[3], bm_sin[3];
	bm_cos[0] = boltzmann_cos(y[0],       THRESHOLD_SLOPE);
	bm_cos[1] = boltzmann_cos(y[N_EQ1],   THRESHOLD_SLOPE);
	bm_cos[2] = boltzmann_cos(y[2*N_EQ1], THRESHOLD_SLOPE);

	bm_sin[0] = boltzmann_sin(y[0],       THRESHOLD_SLOPE);
	bm_sin[1] = boltzmann_sin(y[N_EQ1],   THRESHOLD_SLOPE);
	bm_sin[2] = boltzmann_sin(y[2*N_EQ1], THRESHOLD_SLOPE);

	derivs_one(y, dxdt, p);
	derivs_one(y+N_EQ1, dxdt+N_EQ1, p+NP1);
	derivs_one(y+2*N_EQ1, dxdt+2*N_EQ1, p+2*NP1);

	dxdt[0]       -= (kij[0]*bm_cos[1] + kij[1]*bm_cos[2])*(1.-2.*bm_sin[0]);
	dxdt[N_EQ1]   -= (kij[2]*bm_cos[0] + kij[3]*bm_cos[2])*(1.-2.*bm_sin[1]);
	dxdt[2*N_EQ1] -= (kij[4]*bm_cos[0] + kij[5]*bm_cos[1])*(1.-2.*bm_sin[2]);
}


// kij[0]  :  1 -o  0
// kij[1]  :  2 -o  0
// kij[2]  :  3 -o  0
// kij[3]  :  0 -o  1
// kij[4]  :  2 -o  1
// kij[5]  :  3 -o  1
// kij[6]  :  0 -o  2
// kij[7]  :  1 -o  2
// kij[8]  :  3 -o  2
// kij[9]  :  0 -o  3
// kij[10] :  1 -o  3
// kij[11] :  2 -o  3
void derivs_four(double* y, double* dxdt, const double* p, const double* kij)
{
	double bm_cos[3], bm_sin[3];
	bm_cos[0] = boltzmann_cos(y[0],       THRESHOLD_SLOPE);
	bm_cos[1] = boltzmann_cos(y[N_EQ1],   THRESHOLD_SLOPE);
	bm_cos[2] = boltzmann_cos(y[2*N_EQ1], THRESHOLD_SLOPE);
	bm_cos[3] = boltzmann_cos(y[3*N_EQ1], THRESHOLD_SLOPE);

	bm_sin[0] = boltzmann_sin(y[0],       THRESHOLD_SLOPE);
	bm_sin[1] = boltzmann_sin(y[N_EQ1],   THRESHOLD_SLOPE);
	bm_sin[2] = boltzmann_sin(y[2*N_EQ1], THRESHOLD_SLOPE);
	bm_sin[3] = boltzmann_sin(y[3*N_EQ1], THRESHOLD_SLOPE);

	derivs_one(y, dxdt, p);
	derivs_one(y+N_EQ1, dxdt+N_EQ1, p+NP1);
	derivs_one(y+2*N_EQ1, dxdt+2*N_EQ1, p+2*NP1);
	derivs_one(y+3*N_EQ1, dxdt+3*N_EQ1, p+3*NP1);

	dxdt[0]       -= (kij[0]*bm_cos[1] + kij[1]* bm_cos[2] + kij[2]* bm_cos[3])*(1.-2.*bm_sin[0]);
	dxdt[N_EQ1]   -= (kij[3]*bm_cos[0] + kij[4]* bm_cos[2] + kij[5]* bm_cos[3])*(1.-2.*bm_sin[1]);
	dxdt[2*N_EQ1] -= (kij[6]*bm_cos[0] + kij[7]* bm_cos[1] + kij[8]* bm_cos[3])*(1.-2.*bm_sin[2]);
	dxdt[3*N_EQ1] -= (kij[9]*bm_cos[0] + kij[10]*bm_cos[1] + kij[11]*bm_cos[2])*(1.-2.*bm_sin[3]);
}


void integrate_three_rk4(double* y, const double* params, const double* kij, double* output, const double dt, const unsigned N, const unsigned stride) {
	unsigned i, j, k;
	double dt2, dt6;
	double y1[N_EQ3], y2[N_EQ3], k1[N_EQ3], k2[N_EQ3], k3[N_EQ3], k4[N_EQ3];

	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<3; j++)
		output[j] = y[j*N_EQ1]; 	

	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_three(y, k1, params, kij);

			for(k=0; k<N_EQ3; k++)
				y1[k] = y[k]+k1[k]*dt2; 			

			derivs_three(y1, k2, params, kij);

			for(k=0; k<N_EQ3; k++)
				y2[k] = y[k]+k2[k]*dt2;

			derivs_three(y2, k3, params, kij);

			for(k=0; k<N_EQ3; k++)
				y2[k] = y[k]+k3[k]*dt; 			

			derivs_three(y2, k4, params, kij);

			for(k=0; k<N_EQ3; k++)
			{
				y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
				y[k] = fmod(y[k], 2.*M_PI);
			}
		}
		for(j=0; j<3; j++)
			output[3*i+j] = y[j*N_EQ1]; 					
	}
}


void integrate_four_rk4(double* y, const double* params, const double* kij, double* output, const double dt, const unsigned N, const unsigned stride) {
	unsigned i, j, k;
	double dt2, dt6;
	double y1[N_EQ4], y2[N_EQ4], k1[N_EQ4], k2[N_EQ4], k3[N_EQ4], k4[N_EQ4];

	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<4; j++)
		output[j] = y[j*N_EQ1]; 	

	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_four(y, k1, params, kij);

			for(k=0; k<N_EQ4; k++)
				y1[k] = y[k]+k1[k]*dt2; 			

			derivs_four(y1, k2, params, kij);

			for(k=0; k<N_EQ4; k++)
				y2[k] = y[k]+k2[k]*dt2;

			derivs_four(y2, k3, params, kij);

			for(k=0; k<N_EQ4; k++)
				y2[k] = y[k]+k3[k]*dt; 			

			derivs_four(y2, k4, params, kij);

			for(k=0; k<N_EQ4; k++)
			{
				y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
				y[k] = fmod(y[k], 2.*M_PI);
			}
		}
		for(j=0; j<4; j++)
			output[4*i+j] = y[j*N_EQ1]; 					
	}
}

}
