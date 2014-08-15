
/*
 * FitzHugh Nagumo oscillator according to Schwabedal and Pikovsky, PRE 2010
 */


#include <math.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#define N_EQ1		2

#define COUPLING_THRESHOLD	0. //
#define THRESHOLD_SLOPE		100. //

extern "C" {

double boltzmann(const double V, const double V_0, const double k)
{
	return 1./(1.+exp(-k*(V-V_0)));
}


void derivs_one(double* y, double* dxdt, const double* p)
{
	dxdt[0] = p[5]*(y[0]-y[0]*y[0]*y[0])-y[1]+p[0]; // x' = m (x-x^3)-y+I
	dxdt[1] = p[1]*(boltzmann(y[0], p[2], p[3])-y[1]); // y' = e (Bfun(x, x_0, k_2) - y)
}



void derivs_n(const unsigned num_osci, double* y, double* dxdt, const double* p, const double g_inh)
{
	unsigned i, j;
	double bm_factor[num_osci], bmf_sum;
	for(i=0; i<num_osci; i++)
	{
		bm_factor[i] = boltzmann(y[i*N_EQ1], COUPLING_THRESHOLD, THRESHOLD_SLOPE);
		derivs_one(y+i*N_EQ1, dxdt+i*N_EQ1, p);
	}
	for(i=0; i<num_osci; i++)
	{
		bmf_sum = 0.;
		for(j=0; j<i; j++)
		{
			bmf_sum += bm_factor[j];
		}
		for(j=i+1; j<num_osci; j++)
		{
			bmf_sum += bm_factor[j];
		}
		dxdt[i*N_EQ1] += g_inh*(p[4]-y[i*N_EQ1])*bmf_sum;
	}
}



void step_n_em(const unsigned num_osci, double* y, const double* params, const double g_inh, const double dt, const unsigned stride, double* noise)
{
	unsigned i, j;
	unsigned n_eq = num_osci*N_EQ1;
	double right_hand_side[n_eq];

	for(i=0; i<stride; i++)
	{
		derivs_n(num_osci, y, right_hand_side, params, g_inh);

		for(j=0; j<n_eq; j++)
		{
			y[j] += right_hand_side[j]*dt; 			
		}
		for(j=0; j<num_osci; j++)
		{
			y[N_EQ1*j] += noise[i*num_osci+j];
		}
	}
}


void integrate_n_em(const unsigned num_osci, double* y, const double* params, const double g_inh, double* output, const double dt, const unsigned N, const unsigned stride, const unsigned seed)
{

	gsl_rng *r;
	unsigned i, j, k;
	double sqdt, noise[num_osci*stride];

	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r, seed);

	sqdt_sigma = *sqrt(dt);

	for(i=0; i<num_osci; i++)
	{
		output[i] = y[i*N_EQ1]; 	
	}

	for(i=1; i<N; i++)
	{
		for(k=0; k<num_osci*stride; k++)
		{
			noise[k] = gsl_ran_gaussian(r, sqdt*sigma);
		}

		step_n_em(num_osci, y, params, g_inh, dt, stride, noise);
		for(j=0; j<num_osci; j++)
		{
			output[num_osci*i+j] = y[j*N_EQ1]; 					
		}
	}
}

}
