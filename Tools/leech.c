
/*
 * Three-cell Leech heart interneuron model is implimented 
 * as in the draft  "Bifurcations of bursting polyrhythms
 * in 3-cell central pattern generators" of Jeremy Wojcik,
 * Robert Clewley, and Andrey Shilnikov, prepared for
 * Journal of Neuroscience 2012.
 */


#include <math.h>
#include <stdio.h>


#define N_EQ1	3
#define N_EQ3	3*N_EQ1
#define N_EQ4	4*N_EQ1

#define THRESHOLD_SLOPE	1000.
#define E_syn	-0.0625  // Sajia: same
#define EmNa	-0.0305  // Sajia: 0.0305
#define EhNa	-0.0325  // Sajia: 0.0333
#define V_K2	-0.018   // Sajia: same


void derivs_one(const double* y, double* dxdt, const double* p)
{
	dxdt[0] = ( -p[5]*pow(1./(1.+exp(-150.*(y[0]-EmNa))), 3.)*y[1]*(y[0]-p[2]) - p[4]*y[2]*y[2]*(y[0]-p[1]) - p[3]*(y[0]-p[0]) - p[9])/p[6]; // dV/dt
	dxdt[1] = (1./(1.+exp(500.*(y[0]-EhNa)))-y[1])/p[7];	// dh_Na/dt
	dxdt[2] = (1./(1.+exp(-83.*(y[0]-V_K2+p[10])))-y[2])/p[8];	// dm_K/dt
}


void derivs_two(const double* y, double* dxdt, const double* p, const double* kij)
{
	double coupling = kij[0]*(E_syn-y[0])/(1.+exp(-1000.*(y[3]-p[11])));
	dxdt[0] = ( -p[5]*pow(1./(1.+exp(-150.*(y[0]-EmNa))), 3.)*y[1]*(y[0]-p[2]) - p[4]*y[2]*y[2]*(y[0]-p[1]) - p[3]*(y[0]-p[0]) - p[9] + coupling)/p[6];
	dxdt[1] = (1./(1.+exp(500.*(y[0]-EhNa)))-y[1])/p[7];
	dxdt[2] = (1./(1.+exp(-83.*(y[0]-V_K2+p[10])))-y[2])/p[8];

	coupling = kij[1]*(E_syn-y[3])/(1.+exp(-1000.*(y[0]-p[11])));
	dxdt[3] = ( -p[5]*pow(1./(1.+exp(-150.*(y[3]-EmNa))), 3.)*y[4]*(y[3]-p[2]) - p[4]*y[5]*y[5]*(y[3]-p[1]) - p[3]*(y[3]-p[0]) - p[9] + coupling )/p[6];
	dxdt[4] = (1./(1.+exp(500.*(y[3]-EhNa)))-y[4])/p[7];
	dxdt[5] = (1./(1.+exp(-83.*(y[3]-V_K2+p[10])))-y[5])/p[8];
}


void derivs_three(const double* y, double* dxdt, const double* p, const double* kij)
{
	double coupling = kij[0]*(E_syn-y[0])/(1.+exp(-1000.*(y[3]-p[11]))) + kij[1]*(E_syn-y[0])/(1.+exp(-1000.*(y[6]-p[11]))) + kij[6]*(y[3]-y[0]) + kij[7]*(y[6]-y[0]);
	dxdt[0] = (-p[5]*pow(1./(1.+exp(-150.*(y[0]-EmNa))), 3.)*y[1]*(y[0]-p[2])-p[4]*y[2]*y[2]*(y[0]-p[1])-p[3]*(y[0]-p[0])-p[9]+coupling)/p[6];
	dxdt[1] = (1./(1.+exp(500.*(y[0]-EhNa)))-y[1])/p[7];													
	dxdt[2] = (1./(1.+exp(-83.*(y[0]-V_K2+p[10])))-y[2])/p[8];

	coupling = kij[2]*(E_syn-y[3])/(1.+exp(-1000.*(y[0]-p[11]))) + kij[3]*(E_syn-y[3])/(1.+exp(-1000.*(y[6]-p[11]))) + kij[6]*(y[0]-y[3]) + kij[8]*(y[6]-y[3]);
	dxdt[3] = ( -p[5]*pow(1./(1.+exp(-150.*(y[3]-EmNa))), 3.)*y[4]*(y[3]-p[2]) - p[4]*y[5]*y[5]*(y[3]-p[1]) - p[3]*(y[3]-p[0]) - p[9] + coupling )/p[6];
	dxdt[4] = (1./(1.+exp(500.*(y[3]-EhNa)))-y[4])/p[7];
	dxdt[5] = (1./(1.+exp(-83.*(y[3]-V_K2+p[10])))-y[5])/p[8];

	coupling = kij[4]*(E_syn-y[6])/(1.+exp(-1000.*(y[0]-p[11]))) + kij[5]*(E_syn-y[6])/(1.+exp(-1000.*(y[3]-p[11]))) + kij[7]*(y[0]-y[6]) + kij[8]*(y[3]-y[6]);
	dxdt[6] = (-p[5]*pow(1./(1.+exp(-150.*(y[6]-EmNa))), 3.)*y[7]*(y[6]-p[2])-p[4]*y[8]*y[8]*(y[6]-p[1])-p[3]*(y[6]-p[0])-p[9]+coupling)/p[6];
	dxdt[7] = (1./(1.+exp(500.*(y[6]-EhNa)))-y[7])/p[7];												     
	dxdt[8] = (1./(1.+exp(-83.*(y[6]-V_K2+p[10])))-y[8])/p[8];											     
}

double boltzmann(const double x, const double x_0, const double k)
{
	return 1./(1.+exp(k*(x_0-x)));
}


void derivs_four(const double* y, double* dxdt, const double* p, const double* kij)
{
	int i;
	double bm_factor[4], V[4];

	for(i=0; i<4; i++)
	{
		V[i] = y[i*N_EQ1];
		derivs_one(y+i*N_EQ1, dxdt+i*N_EQ1, p);
		bm_factor[i] = boltzmann(V[i], p[11], THRESHOLD_SLOPE)/p[6];
	}

	dxdt[0] +=       (E_syn-V[0])*(kij[0]*bm_factor[1] + kij[1]* bm_factor[2] + kij[2]* bm_factor[3]);
	dxdt[N_EQ1] +=   (E_syn-V[1])*(kij[3]*bm_factor[0] + kij[4]* bm_factor[2] + kij[5]* bm_factor[3]);
	dxdt[2*N_EQ1] += (E_syn-V[2])*(kij[6]*bm_factor[0] + kij[7]* bm_factor[1] + kij[8]* bm_factor[3]);
	dxdt[3*N_EQ1] += (E_syn-V[3])*(kij[9]*bm_factor[0] + kij[10]*bm_factor[1] + kij[11]*bm_factor[2]);

	dxdt[0]       += (kij[12]*(V[1]-V[0]) + kij[13]*(V[2]-V[0]) + kij[14]*(V[3]-V[0]))/p[6];
	dxdt[N_EQ1]   += (kij[12]*(V[0]-V[1]) + kij[15]*(V[2]-V[1]) + kij[16]*(V[3]-V[1]))/p[6];
	dxdt[2*N_EQ1] += (kij[13]*(V[0]-V[2]) + kij[15]*(V[1]-V[2]) + kij[17]*(V[3]-V[2]))/p[6];
	dxdt[3*N_EQ1] += (kij[14]*(V[0]-V[3]) + kij[16]*(V[1]-V[3]) + kij[17]*(V[2]-V[3]))/p[6];
}


void integrate_four_rk4(double* y, const double* params, const double* coupling, double* output, const double dt, const unsigned N, const unsigned stride)
{
	unsigned i, j, k;
	double dt2, dt6;
	double y1[N_EQ4], y2[N_EQ4], k1[N_EQ4], k2[N_EQ4], k3[N_EQ4], k4[N_EQ4];
	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<4; j++)
		output[j] = y[3*j];

	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_four(y, k1, params, coupling);
			for(k=0; k<N_EQ4; k++)
				y1[k] = y[k]+k1[k]*dt2; 			

			derivs_four(y1, k2, params, coupling);

			for(k=0; k<N_EQ4; k++)
				y2[k] = y[k]+k2[k]*dt2; 			
			derivs_four(y2, k3, params, coupling);

			for(k=0; k<N_EQ4; k++)
				y2[k] = y[k]+k3[k]*dt; 			

			derivs_four(y2, k4, params, coupling);
			for(k=0; k<N_EQ4; k++)
				y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
		}
		for(j=0; j<4; j++)
			output[4*i+j] = y[3*j]; 					
	}
};


void integrate_three_rk4(double* y, const double* params, const double* coupling, double* output, const double dt, const unsigned N, const unsigned stride)
{
	unsigned i, j, k;
	double dt2, dt6;
	double y1[N_EQ3], y2[N_EQ3], k1[N_EQ3], k2[N_EQ3], k3[N_EQ3], k4[N_EQ3];
	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<3; j++)
		output[j] = y[3*j]; 	

	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_three(y, k1, params, coupling);
			for(k=0; k<N_EQ3; k++)
				y1[k] = y[k]+k1[k]*dt2; 			

			derivs_three(y1, k2, params, coupling);

			for(k=0; k<N_EQ3; k++)
				y2[k] = y[k]+k2[k]*dt2; 			
			derivs_three(y2, k3, params, coupling);

			for(k=0; k<N_EQ3; k++)
				y2[k] = y[k]+k3[k]*dt; 			

			derivs_three(y2, k4, params, coupling);
			for(k=0; k<N_EQ3; k++)
				y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
		}
		for(j=0; j<3; j++)
			output[3*i+j] = y[3*j]; 					
	}
};


void integrate_one_rk4(double* y, const double dt, const unsigned N, const unsigned stride, const double* P, double* output)
{
	unsigned i, j, k;
	double dt2, dt6;
	double y1[3], y2[3], k1[3], k2[3], k3[3], k4[3];
	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<N_EQ1; j++)
		output[j] = y[j];



	for(i=1; i<N; i++)
	{
		for(j=0; j<stride; j++)
		{
			derivs_one(y, k1, P);
			for(k=0; k<N_EQ1; k++)
			       y1[k] = y[k]+k1[k]*dt2; 			

			derivs_one(y1, k2, P);
			for(k=0; k<N_EQ1; k++)
				y2[k] = y[k]+k2[k]*dt2; 			

			derivs_one(y2, k3, P);
			for(k=0; k<N_EQ1; k++)
				y2[k] = y[k]+k3[k]*dt; 			

			derivs_one(y2, k4, P);
			for(k=0; k<N_EQ1; k++)
				y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
		}
		for(j=0; j<N_EQ1; j++) output[N_EQ1*i+j] = y[j];
	}
}



void step_three_rk4(double* y, const double* params, const double* coupling, const double dt, const unsigned stride)
{
	unsigned i, j, k;
	double dt2, dt6;
	double y1[N_EQ3], y2[N_EQ3], k1[N_EQ3], k2[N_EQ3], k3[N_EQ3], k4[N_EQ3];
	dt2 = dt/2.; dt6 = dt/6.;

	for(j=0; j<stride; j++)
	{
		derivs_three(y, k1, params, coupling);
		for(k=0; k<N_EQ3; k++)
			y1[k] = y[k]+k1[k]*dt2; 			

		derivs_three(y1, k2, params, coupling);

		for(k=0; k<N_EQ3; k++)
			y2[k] = y[k]+k2[k]*dt2; 			
		derivs_three(y2, k3, params, coupling);

		for(k=0; k<N_EQ3; k++)
			y2[k] = y[k]+k3[k]*dt; 			

		derivs_three(y2, k4, params, coupling);
		for(k=0; k<N_EQ3; k++)
			y[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
	}
};








