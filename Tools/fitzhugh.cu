
#include<iostream>
#include<cstdlib>
#include<vector>
#include<cuda.h>

using namespace std;

#define THREADS_PER_BLOCK 32 // threads per block

#define N_EQ1		2
#define NP1		6

#define MAX_ITERATIONS	10000000 // 10^8

#define NODES3		3
#define N_EQ3		(NODES3*N_EQ1)
#define NP3		(NODES3*NP1)
#define NK3		((NODES3-1)*NODES3)

#define NODES4		4
#define N_EQ4		(NODES4*N_EQ1)
#define NP4		(NODES4*NP1)
#define NK4		((NODES4-1)*NODES4)

#define COUPLING_THRESHOLD	0. //
#define THRESHOLD_SLOPE		100. //


__device__ double boltzmann(const double V, const double V_0, const double k)
{
	return 1./(1.+exp(-k*(V-V_0)));
}

__device__ void derivs_one(double *y, double *dxdt, const double *p)
{
	dxdt[0] = p[5]*(y[0]-y[0]*y[0]*y[0])-y[1]+p[0]; // x' = m (x-x^3)-y+I
	dxdt[1] = p[1]*(boltzmann(y[0], p[2], p[3])-y[1]); // y' = e (Bfun(x, x_0, k_2) - y)
}


__device__ void derivs_three(double* y, double* dxdt, const double* p, const double* kij)
{
	double bm_factor[3];
	bm_factor[0] = boltzmann(y[0],       COUPLING_THRESHOLD, THRESHOLD_SLOPE);
	bm_factor[1] = boltzmann(y[N_EQ1],   COUPLING_THRESHOLD, THRESHOLD_SLOPE);
	bm_factor[2] = boltzmann(y[2*N_EQ1], COUPLING_THRESHOLD, THRESHOLD_SLOPE);

	derivs_one(y, dxdt, p);
	derivs_one(y+N_EQ1, dxdt+N_EQ1, p+NP1);
	derivs_one(y+2*N_EQ1, dxdt+2*N_EQ1, p+2*NP1);

	dxdt[0]       += (p[4]-y[0])*            (kij[0]*bm_factor[1] + kij[1]*bm_factor[2]);
	dxdt[N_EQ1]   += (p[4+NP1]-y[N_EQ1])*    (kij[2]*bm_factor[0] + kij[3]*bm_factor[2]);
	dxdt[2*N_EQ1] += (p[4+2*NP1]-y[2*N_EQ1])*(kij[4]*bm_factor[0] + kij[5]*bm_factor[1]);
}
 

__device__ void assign(double *x, double *y, const unsigned ngls) { for(int i=0; i<ngls; i++) x[i] = y[i]; }


__global__ void rk4_three(const unsigned THREADS, double *X, const double dt, const unsigned iterations, const unsigned stride, double *p, double *kij, double *output)
{
	const int THREAD_ID=blockDim.x*blockIdx.x + threadIdx.x;

	if(THREAD_ID < THREADS)
	{

	const int thread_init=N_EQ3*THREAD_ID;
	const double dt2=dt/2., dt6=dt/6.;
	double x[N_EQ3], x1[N_EQ3], x2[N_EQ3], k1[N_EQ3], k2[N_EQ3], k3[N_EQ3], k4[N_EQ3], P[NP3], K[NK3];
	int t, s, i;

	assign(x, X+thread_init, N_EQ3);
	assign(P, p, NP3);
	assign(K, kij, NK3);

	for(i=0; i<NODES3; i++)
		output[THREAD_ID*NODES3*iterations + i] = x[N_EQ1*i];	// save Voltages (x_0, x_3, x_6) to output[THREAD_ID*3:(THREAD_ID+1)*3]

	for(t=1; t<iterations; t++)
	{
		for(s=0; s<stride; s++)
		{
			derivs_three(x, k1, P, K); 							// k1 = f(x)
			for(i=0; i<N_EQ3; i++) x1[i] = x[i]+k1[i]*dt2; 				// x1 = x + k1*dt/2
			derivs_three(x1, k2, P, K); 							// k2 = f(x1)
			for(i=0; i<N_EQ3; i++) x2[i] = x[i]+k2[i]*dt2; 				// x2 = x + k2*dt/2
			derivs_three(x2, k3, P, K); 							// k3 = f(x2)
			for(i=0; i<N_EQ3; i++) x2[i] = x[i]+k3[i]*dt; 				// x2 = x + k3*dt
			derivs_three(x2, k4, P, K); 							// k4 = f(x+k3*dt)
			for(i=0; i<N_EQ3; i++) x[i] += dt6*(k1[i]+2.*(k2[i]+k3[i])+k4[i]); 	// x_n+1 = x_n + dt (...)/6.
		}
		for(s=0; s<NODES3; s++)
			output[THREAD_ID*NODES3*iterations + NODES3*t + s] = x[N_EQ1*s];	// save Voltages (x_0, x_3, x_6) to output[THREAD_ID*3:(THREAD_ID+1)*3]
	}

	}
};


__device__ void rk4_three_step(double *x, const double *P, const double *K, const double dt, const unsigned stride)
{
	const double dt2=0.5*dt, dt6=dt/6.;
	double x1[N_EQ3], x2[N_EQ3], k1[N_EQ3], k2[N_EQ3], k3[N_EQ3], k4[N_EQ3];
	int s, i;
	for(s=0; s<stride; s++)
	{
		derivs_three(x, k1, P, K); 							// k1 = f(x)
		for(i=0; i<N_EQ3; i++)
		{
			x1[i] = x[i]+k1[i]*dt2; 						// x1 = x + k1*dt/2
		}
		derivs_three(x1, k2, P, K); 							// k2 = f(x1)
		for(i=0; i<N_EQ3; i++)
		{
			x2[i] = x[i]+k2[i]*dt2; 						// x2 = x + k2*dt/2
		}
		derivs_three(x2, k3, P, K); 							// k3 = f(x2)
		for(i=0; i<N_EQ3; i++)
		{
			x2[i] = x[i]+k3[i]*dt; 							// x2 = x + k3*dt
		}
		derivs_three(x2, k4, P, K); 							// k4 = f(x+k3*dt)
		for(i=0; i<N_EQ3; i++)
		{
			x[i] += dt6*(k1[i]+2.*(k2[i]+k3[i])+k4[i]); 				// x_n+1 = x_n + dt (...)/6.
		}
	}
}



__global__ void rk4_three_save_crossings(const unsigned THREADS,
					double *X, const double dt,
					const double threshold, const unsigned num_crossings,
					const unsigned stride, double *p, double *kij, double *output)
{
	const int THREAD_ID=blockDim.x*blockIdx.x + threadIdx.x;

	if(THREAD_ID < THREADS)
	{

	const int thread_init=N_EQ3*THREAD_ID;
	const double DT = dt*(double)stride;
	double current, previous[NODES3], x[N_EQ3], P[NP3], K[NK3];
	int t, i, i_c[NODES3];
	bool done, condition;

	assign(x, X+thread_init, N_EQ3);
	assign(P, p, NP3);
	assign(K, kij, NK3);

	for(i=0; i<NODES3; i++)
		i_c[i] = 0;

	for(t=1; t<MAX_ITERATIONS; t++)
	{
		rk4_three_step(x, P, K, dt, stride);

		done = true;

		for(i=0; i<NODES3; i++)
		{
			condition = i_c[i] >= num_crossings;
			done = done && condition; 	// false if not done yet.

			if(!condition)
			{
				current = x[N_EQ1*i];
				output[THREAD_ID*NODES3*num_crossings + NODES3*i_c[i] + i] = (double)t*DT;
				i_c[i] += (int)( previous[i]<threshold && current>threshold );
				previous[i] = current;
			}
		}

		if(done) break;
	}

	}
};



__device__ void derivs_four(double* y, double* dxdt, const double* p, const double* kij)
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


__global__ void rk4_four(const unsigned THREADS, double *X, const double dt, const unsigned iterations, const unsigned stride, double *p, double *kij, double *output)
{
	const int THREAD_ID=blockDim.x*blockIdx.x + threadIdx.x;

	if(THREAD_ID < THREADS)
	{

	const int thread_init=N_EQ4*THREAD_ID;
	const double dt2=dt/2., dt6=dt/6.;
	double x[N_EQ4], x1[N_EQ4], x2[N_EQ4], k1[N_EQ4], k2[N_EQ4], k3[N_EQ4], k4[N_EQ4], P[NP4], K[NK4];
	int t, s, i;

	assign(x, X+thread_init, N_EQ4);
	assign(P, p, NP4);
	assign(K, kij, NK4);

	for(i=0; i<NODES4; i++)
		output[THREAD_ID*NODES4*iterations + i] = x[N_EQ1*i];	// save Voltages (x_0, x_3, x_6) to output[THREAD_ID*3:(THREAD_ID+1)*3]

	for(t=1; t<iterations; t++)
	{
		for(s=0; s<stride; s++)
		{
			derivs_four(x, k1, P, K); 							// k1 = f(x)
			for(i=0; i<N_EQ4; i++) x1[i] = x[i]+k1[i]*dt2; 				// x1 = x + k1*dt/2
			derivs_four(x1, k2, P, K); 							// k2 = f(x1)
			for(i=0; i<N_EQ4; i++) x2[i] = x[i]+k2[i]*dt2; 				// x2 = x + k2*dt/2
			derivs_four(x2, k3, P, K); 							// k3 = f(x2)
			for(i=0; i<N_EQ4; i++) x2[i] = x[i]+k3[i]*dt; 				// x2 = x + k3*dt
			derivs_four(x2, k4, P, K); 							// k4 = f(x+k3*dt)
			for(i=0; i<N_EQ4; i++) x[i] += dt6*(k1[i]+2.*(k2[i]+k3[i])+k4[i]); 	// x_n+1 = x_n + dt (...)/6.
		}
		for(s=0; s<NODES4; s++)
			output[THREAD_ID*NODES4*iterations + NODES4*t + s] = x[N_EQ1*s];	// save Voltages (x_0, x_3, x_6) to output[THREAD_ID*3:(THREAD_ID+1)*3]
	}

	}
};


extern "C" {


void cuda_integrate_three(double* y_init, const unsigned size_y_init,
		     double* y_out,
		     double* coupling,
		     double* params,
		const double dt, const unsigned iterations, const unsigned stride)
{
	unsigned size_y_out, THREADS, BLOCKS;
	THREADS = size_y_init/(NODES3*N_EQ1);
	BLOCKS = THREADS/THREADS_PER_BLOCK+1;
	size_y_out = iterations*THREADS*NODES3;


	cudaError_t err;
	int device=0;
	cudaDeviceProp device_prop;

	if(err = cudaGetDeviceProperties(&device_prop, device)) cout << cudaGetErrorString(err) << endl;
        cout << "# GPU Device " << device << ": \""<< device_prop.name << "\" with compute capability " << device_prop.major << "." << device_prop.minor << endl;

	double *p, *kij, *Xd, *output;

	if(err = cudaMalloc((void**) &p, NP3*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &kij, NK3*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &Xd, size_y_init*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &output, size_y_out*sizeof(double))) cout << cudaGetErrorString(err) << endl;



	if(err = cudaMemcpy(p, params, NP3*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMemcpy(kij, coupling, NK3*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMemcpy(Xd, y_init, size_y_init*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;


	cout << "# starting computation with " << BLOCKS << " blocks a " << THREADS_PER_BLOCK << " threads." << endl;
	rk4_three<<<BLOCKS, THREADS_PER_BLOCK>>>(THREADS, Xd, dt, iterations, stride, p, kij, output);


	if(err = cudaMemcpy(y_out, output, (size_t)size_y_out*sizeof(double), cudaMemcpyDeviceToHost)) cout << "Error in copy device to host: " << cudaGetErrorString(err) << endl;


	cudaFree(p);
	cudaFree(kij);
	cudaFree(Xd);
	cudaFree(output);
}


void cuda_crossing_three(double* y_init, const unsigned size_y_init,
		     double* y_out,
		     double* coupling,
		     double* params,
		const double threshold, const unsigned num_crossings,
	       	const double dt, const unsigned stride)
{
	unsigned size_y_out, THREADS, BLOCKS;
	THREADS = size_y_init/N_EQ3;
	BLOCKS = THREADS/THREADS_PER_BLOCK+1;
	size_y_out = num_crossings*THREADS*NODES3;


	cudaError_t err;
	int device=0;
	cudaDeviceProp device_prop;

	if(err = cudaGetDeviceProperties(&device_prop, device)) cout << cudaGetErrorString(err) << endl;
        cout << "# GPU Device " << device << ": \""<< device_prop.name << "\" with compute capability " << device_prop.major << "." << device_prop.minor << endl;


	double *p, *kij, *Xd, *output;

	if(err = cudaMalloc((void**) &p, NP4*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &kij, NK4*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &Xd, size_y_init*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &output, size_y_out*sizeof(double))) cout << cudaGetErrorString(err) << endl;

	if(err = cudaMemcpy(p, params, NP4*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMemcpy(kij, coupling, NK4*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMemcpy(Xd, y_init, size_y_init*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;


	cout << "# starting computation with " << BLOCKS << " blocks a " << THREADS_PER_BLOCK << " threads." << endl;

	rk4_three_save_crossings<<<BLOCKS, THREADS_PER_BLOCK>>>(THREADS, Xd, dt, threshold, num_crossings, stride, p, kij, output);

	if(err = cudaMemcpy(y_out, output, (size_t)size_y_out*sizeof(double), cudaMemcpyDeviceToHost)) cout << "Error in copy device to host: " << cudaGetErrorString(err) << endl;

	cudaFree(p);
	cudaFree(kij);
	cudaFree(Xd);
	cudaFree(output);
}


void cuda_integrate_four(double* y_init, const unsigned size_y_init,
		     double* y_out,
		     double* coupling,
		     double* params,
		const double dt, const unsigned iterations, const unsigned stride)
{
	unsigned size_y_out, THREADS, BLOCKS;
	THREADS = size_y_init/N_EQ4;
	BLOCKS = THREADS/THREADS_PER_BLOCK+1;
	size_y_out = iterations*THREADS*NODES4;


	cudaError_t err;
	int device=0;
	cudaDeviceProp device_prop;

	if(err = cudaGetDeviceProperties(&device_prop, device)) cout << cudaGetErrorString(err) << endl;
        cout << "# GPU Device " << device << ": \""<< device_prop.name << "\" with compute capability " << device_prop.major << "." << device_prop.minor << endl;


	double *p, *kij, *Xd, *output;

	if(err = cudaMalloc((void**) &p, NP4*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &kij, NK4*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &Xd, size_y_init*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &output, size_y_out*sizeof(double))) cout << cudaGetErrorString(err) << endl;

	if(err = cudaMemcpy(p, params, NP4*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMemcpy(kij, coupling, NK4*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMemcpy(Xd, y_init, size_y_init*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;




	cout << "# starting computation with " << BLOCKS << " blocks a " << THREADS_PER_BLOCK << " threads." << endl;
	rk4_four<<<BLOCKS, THREADS_PER_BLOCK>>>(THREADS, Xd, dt, iterations, stride, p, kij, output);


	if(err = cudaMemcpy(y_out, output, (size_t)size_y_out*sizeof(double), cudaMemcpyDeviceToHost)) cout << "Error in copy device to host: " << cudaGetErrorString(err) << endl;


	cudaFree(p);
	cudaFree(kij);
	cudaFree(Xd);
	cudaFree(output);
}


}
