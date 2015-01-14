
/*
 * FitzHugh Nagumo oscillator according to Schwabedal and Pikovsky, PRE 2010
 */


#include <math.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#define N_EQ1		2
#define NP1		6
#define N_EQ3		(3*N_EQ1)
#define NP3		(3*NP1)	// system parameters + noise strength
#define N_EQ4		(4*N_EQ1)
#define NP4		(4*NP1)
#define NC		6


#define COUPLING_THRESHOLD	0. //
#define THRESHOLD_SLOPE		100. //

double boltzmann(const double V, const double V_0, const double k)
{
	return 1./(1.+exp(-k*(V-V_0)));
}


void derivs_one(double* y, double* dxdt, const double* p)
{
	dxdt[0] = p[5]*(y[0]-y[0]*y[0]*y[0])-y[1]+p[0]; // x' = m (x-x^3)-y+I
	dxdt[1] = p[1]*(boltzmann(y[0], p[2], p[3])-y[1]); // y' = e (Bfun(x, x_0, k_2) - y)
}


void integrate_one_em(double* y, const double* params, double* output, const double dt, const unsigned N, const unsigned stride, const unsigned seed)
{
	gsl_rng *r;
	unsigned i, j, k;
	double sigma, right_hand_side[N_EQ1];

	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r, seed);

	sigma = params[NP1]*sqrt(dt);	// should be sigma * dt**(1/2)
	for(j=0; j<N_EQ1; j++)
	{
		output[j] = y[j];
	}
	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_one(y, right_hand_side, params);
			for(k=0; k<N_EQ1; k++)
			{
				y[k] += right_hand_side[k]*dt;
			}
			y[0] += gsl_ran_gaussian(r, sigma);
		}
		for(j=0; j<N_EQ1; j++)
		{
			output[N_EQ1*i+j] = y[j];
		}
	}
	gsl_rng_free(r);
}


void derivs_three(double* y, double* dxdt, const double* p, const double* kij)
{
	double bm_factor[3];
	bm_factor[0] = boltzmann(y[0],       COUPLING_THRESHOLD, THRESHOLD_SLOPE);
	bm_factor[1] = boltzmann(y[N_EQ1],   COUPLING_THRESHOLD, THRESHOLD_SLOPE);
	bm_factor[2] = boltzmann(y[2*N_EQ1], COUPLING_THRESHOLD, THRESHOLD_SLOPE);

	derivs_one(y, dxdt, p);
	derivs_one(y+N_EQ1, dxdt+N_EQ1, p+NP1);
	derivs_one(y+2*N_EQ1, dxdt+2*N_EQ1, p+2*NP1);

	dxdt[0] +=       (p[4]-y[0])*            (kij[0]*bm_factor[1] + kij[1]*bm_factor[2]);
	dxdt[N_EQ1] +=   (p[4+NP1]-y[N_EQ1])*    (kij[2]*bm_factor[0] + kij[3]*bm_factor[2]);
	dxdt[2*N_EQ1] += (p[4+2*NP1]-y[2*N_EQ1])*(kij[4]*bm_factor[0] + kij[5]*bm_factor[1]);
}


void derivs_four(double* y, double* dxdt, const double* p, const double* kij)
{
	double bm_factor[4];
	bm_factor[0] = boltzmann(y[0],       COUPLING_THRESHOLD, THRESHOLD_SLOPE);
	bm_factor[1] = boltzmann(y[N_EQ1],   COUPLING_THRESHOLD, THRESHOLD_SLOPE);
	bm_factor[2] = boltzmann(y[2*N_EQ1], COUPLING_THRESHOLD, THRESHOLD_SLOPE);
	bm_factor[3] = boltzmann(y[3*N_EQ1], COUPLING_THRESHOLD, THRESHOLD_SLOPE);

	derivs_one(y, dxdt, p);
	derivs_one(y+N_EQ1, dxdt+N_EQ1, p+NP1);
	derivs_one(y+2*N_EQ1, dxdt+2*N_EQ1, p+2*NP1);
	derivs_one(y+3*N_EQ1, dxdt+3*N_EQ1, p+3*NP1);

	dxdt[0] +=       (p[4]-y[0])*            (kij[0]*bm_factor[1] + kij[1]* bm_factor[2] + kij[2]* bm_factor[3]);
	dxdt[N_EQ1] +=   (p[4+NP1]-y[N_EQ1])*    (kij[3]*bm_factor[0] + kij[4]* bm_factor[2] + kij[5]* bm_factor[3]);
	dxdt[2*N_EQ1] += (p[4+2*NP1]-y[2*N_EQ1])*(kij[6]*bm_factor[0] + kij[7]* bm_factor[1] + kij[8]* bm_factor[3]);
	dxdt[3*N_EQ1] += (p[4+3*NP1]-y[3*N_EQ1])*(kij[9]*bm_factor[0] + kij[10]*bm_factor[1] + kij[11]*bm_factor[2]);
}


void integrate_three_em(double* y, const double* params, const double* kij, double* output, const double dt, const unsigned N, const unsigned stride, const unsigned seed)
{
	gsl_rng *r;
	unsigned i, j, k;
	double sigma, right_hand_side[N_EQ3];

	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r, seed);

	sigma = params[NP3]*sqrt(dt);	// should be sigma * dt**(1/2)

	for(j=0; j<3; j++)
	{
		output[j] = y[j*N_EQ1];
	}
	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_three(y, right_hand_side, params, kij);

			for(k=0; k<N_EQ3; k++)
			{
				y[k] += right_hand_side[k]*dt;	
			}
			for(k=0; k<3; k++)
			{
				y[k*N_EQ1] += gsl_ran_gaussian(r, sigma);
			}
		}
		for(j=0; j<3; j++)
		{
			output[3*i+j] = y[j*N_EQ1]; 					
		}
	}
	gsl_rng_free(r);
}


void integrate_four_em(double* y, const double* params, const double* kij, double* output, const double dt, const unsigned N, const unsigned stride, const unsigned seed)
{
	gsl_rng *r;
	unsigned i, j, k;
	double sigma, right_hand_side[N_EQ4];

	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r, seed);

	sigma = params[NP4]*sqrt(dt);	// should be sigma * dt**(1/2)

	for(j=0; j<4; j++)
	{
		output[j] = y[j*N_EQ1]; 	
	}

	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_four(y, right_hand_side, params, kij);

			for(k=0; k<N_EQ4; k++)
			{
				y[k] += right_hand_side[k]*dt; 			
			}
			for(k=0; k<4; k++)
			{
				y[k*N_EQ1] += gsl_ran_gaussian(r, sigma);
			}
		}
		for(j=0; j<4; j++)
		{
			output[4*i+j] = y[j*N_EQ1]; 					
		}
	}
	gsl_rng_free(r);
}




