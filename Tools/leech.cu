
#include<iostream>
#include<cstdlib>
#include<vector>
#include<cuda.h>

using namespace std;

#define THREADS_PER_BLOCK 32 // threads per block

#define NGLS1	3 	// number of equations in one reduced leech models
#define NP1	11 	// number of parameters per reduced leech model

#define NGLS	9 	// number of equations in 3 reduced leech models
#define NP3	12 	// number of parameters per reduced leech model
#define NK	9 	// number of coupling konstants for reduced leech network


#define Esyn	-0.0625
#define EmNa	-0.0305 
#define EhNa	-0.0325 
#define V_K2	-0.018   // Sajia: same


__device__ void leech_one(double *y, double *dxdt, double *p) {
	dxdt[0] = (-p[5]*pow(1./(1.+exp(-150.*(y[0]-EmNa))), 3)*y[1]*(y[0]-p[2])-p[4]*y[2]*y[2]*(y[0]-p[1])-p[3]*(y[0]-p[0])-p[9])/p[6];
	dxdt[1] = (1./(1.+exp(500.*(y[0]-EhNa)))-y[1])/p[7];
	dxdt[2] = (1./(1.+exp(-83.*(y[0]-V_K2+p[10])))-y[2])/p[8];
};


__device__ void derivs_three(const double* y, double* dxdt, const double* p, const double* kij)
{
	double coupling = kij[0]*(Esyn-y[0])/(1.+exp(-1000.*(y[3]-p[11]))) + kij[1]*(Esyn-y[0])/(1.+exp(-1000.*(y[6]-p[11]))) + kij[6]*(y[3]-y[0]) + kij[7]*(y[6]-y[0]);
	dxdt[0] = (-p[5]*pow(1./(1.+exp(-150.*(y[0]-EmNa))), 3.)*y[1]*(y[0]-p[2])-p[4]*y[2]*y[2]*(y[0]-p[1])-p[3]*(y[0]-p[0])-p[9]+coupling)/p[6];
	dxdt[1] = (1./(1.+exp(500.*(y[0]-EhNa)))-y[1])/p[7];													
	dxdt[2] = (1./(1.+exp(-83.*(y[0]-V_K2+p[10])))-y[2])/p[8];

	coupling = kij[2]*(Esyn-y[3])/(1.+exp(-1000.*(y[0]-p[11]))) + kij[3]*(Esyn-y[3])/(1.+exp(-1000.*(y[6]-p[11]))) + kij[6]*(y[0]-y[3]) + kij[8]*(y[6]-y[3]);
	dxdt[3] = ( -p[5]*pow(1./(1.+exp(-150.*(y[3]-EmNa))), 3.)*y[4]*(y[3]-p[2]) - p[4]*y[5]*y[5]*(y[3]-p[1]) - p[3]*(y[3]-p[0]) - p[9] + coupling )/p[6];
	dxdt[4] = (1./(1.+exp(500.*(y[3]-EhNa)))-y[4])/p[7];
	dxdt[5] = (1./(1.+exp(-83.*(y[3]-V_K2+p[10])))-y[5])/p[8];

	coupling = kij[4]*(Esyn-y[6])/(1.+exp(-1000.*(y[0]-p[11]))) + kij[5]*(Esyn-y[6])/(1.+exp(-1000.*(y[3]-p[11]))) + kij[7]*(y[0]-y[6]) + kij[8]*(y[3]-y[6]);
	dxdt[6] = (-p[5]*pow(1./(1.+exp(-150.*(y[6]-EmNa))), 3.)*y[7]*(y[6]-p[2])-p[4]*y[8]*y[8]*(y[6]-p[1])-p[3]*(y[6]-p[0])-p[9]+coupling)/p[6];	     
	dxdt[7] = (1./(1.+exp(500.*(y[6]-EhNa)))-y[7])/p[7];												     
	dxdt[8] = (1./(1.+exp(-83.*(y[6]-V_K2+p[10])))-y[8])/p[8];											     
}



__device__ void assign(double *x, double *y, const unsigned ngls)
{
	for(unsigned i=0; i<ngls; i++)
	{
		x[i] = y[i];
	}
}



__global__ void rk4_one_relax(const unsigned THREADS, double *X, const double dt, const unsigned ITERATIONS, double *p)
{
	const int THREAD_ID=blockDim.x*blockIdx.x + threadIdx.x;

	if(THREAD_ID < THREADS)
	{

	const int thread_idx=NGLS1*THREAD_ID;
	const double dt2=dt/2., dt6=dt/6.;
	double x[NGLS1], x1[NGLS1], x2[NGLS1], k1[NGLS1], k2[NGLS1], k3[NGLS1], k4[NGLS1], P[NP1];
	int t, i;

	assign(x, X+thread_idx, NGLS1);
	assign(P, p, NP1);

	for(t=0; t<ITERATIONS; t++)
	{
		leech_one(x, k1, P); 							// k1 = f(x)
		for(i=0; i<NGLS1; i++) x1[i] = x[i]+k1[i]*dt2; 				// x1 = x + k1*dt/2
		leech_one(x1, k2, P); 							// k2 = f(x1)
		for(i=0; i<NGLS1; i++) x2[i] = x[i]+k2[i]*dt2; 				// x2 = x + k2*dt/2
		leech_one(x2, k3, P); 							// k3 = f(x2)
		for(i=0; i<NGLS1; i++) x2[i] = x[i]+k3[i]*dt; 				// x2 = x + k3*dt
		leech_one(x2, k4, P); 							// k4 = f(x+k3*dt)
		for(i=0; i<NGLS1; i++) x[i] += dt6*(k1[i]+2.*(k2[i]+k3[i])+k4[i]); 	// x_n+1 = x_n + dt (...)/6.
	}
	assign(X+thread_idx, x, NGLS1);

	}
};


__global__ void rk4_three(const unsigned THREADS, double *X, const double dt, const unsigned ITERATIONS, const unsigned stride, double *p, double *kij, double *output)
{
	const int THREAD_ID=blockDim.x*blockIdx.x + threadIdx.x;
	if(THREAD_ID < THREADS)
	{

	const int thread_idx=NGLS*THREAD_ID;
	const double dt2=dt/2., dt6=dt/6.;
	double x[NGLS], x1[NGLS], x2[NGLS], k1[NGLS], k2[NGLS], k3[NGLS], k4[NGLS], P[NP3], K[NK];
	int t, s, i;

	assign(x, X+thread_idx, NGLS);
	assign(P, p, NP3);
	assign(K, kij, NK);

	for(i=0; i<3; i++)
		output[THREAD_ID*3*ITERATIONS + i] = x[3*i];	// save Voltages (x_0, x_3, x_6) to output[THREAD_ID*3:(THREAD_ID+1)*3]

	for(t=1; t<ITERATIONS; t++)
	{
		for(s=0; s<stride; s++)
		{
			derivs_three(x, k1, P, K); 							// k1 = f(x)
			for(i=0; i<NGLS; i++) x1[i] = x[i]+k1[i]*dt2; 				// x1 = x + k1*dt/2
			derivs_three(x1, k2, P, K); 							// k2 = f(x1)
			for(i=0; i<NGLS; i++) x2[i] = x[i]+k2[i]*dt2; 				// x2 = x + k2*dt/2
			derivs_three(x2, k3, P, K); 							// k3 = f(x2)
			for(i=0; i<NGLS; i++) x2[i] = x[i]+k3[i]*dt; 				// x2 = x + k3*dt
			derivs_three(x2, k4, P, K); 							// k4 = f(x+k3*dt)
			for(i=0; i<NGLS; i++) x[i] += dt6*(k1[i]+2.*(k2[i]+k3[i])+k4[i]); 	// x_n+1 = x_n + dt (...)/6.
		}
		for(s=0; s<3; s++)
			output[THREAD_ID*3*ITERATIONS + 3*t + s] = x[3*s];	// save Voltages (x_0, x_3, x_6) to output[THREAD_ID*3:(THREAD_ID+1)*3]
	}

	}
};



extern "C" {


void cuda_relax_one(double* y_init, const unsigned size_y_init,
		     double* params,
			const double dt, const unsigned ITERATIONS)
{
	unsigned THREADS, BLOCKS;
	THREADS = size_y_init/NGLS1;
	BLOCKS = THREADS/THREADS_PER_BLOCK+1;

	cudaError_t err;
	int device=0;
	cudaDeviceProp device_prop;


	if(err = cudaGetDeviceProperties(&device_prop, device)) cout << cudaGetErrorString(err) << endl;
        cout << "# GPU Device " << device << ": \""<< device_prop.name << "\" with compute capability " << device_prop.major << "." << device_prop.minor << endl;


	double *p, *Xd;
	if(err = cudaMalloc((void**) &p, NP1*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &Xd, size_y_init*sizeof(double))) cout << cudaGetErrorString(err) << endl;

	if(err = cudaMemcpy(p, params, NP1*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMemcpy(Xd, y_init, size_y_init*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;


	rk4_one_relax<<<BLOCKS, THREADS_PER_BLOCK>>>(THREADS, Xd, dt, ITERATIONS, p);


	if(err = cudaMemcpy(y_init, Xd, (size_t)size_y_init*sizeof(double), cudaMemcpyDeviceToHost)) cout << cudaGetErrorString(err) << endl;
	
	cudaFree(p);
	cudaFree(Xd);
}



void cuda_integrate_three(double* y_init, const unsigned size_y_init,
		     double* y_out,
		     double* coupling,
		     double* params,
		const double dt, const unsigned ITERATIONS, const unsigned stride)
{
	unsigned size_y_out, THREADS, BLOCKS, NODES=3;
	THREADS = size_y_init/(NODES*NGLS1);
	BLOCKS = THREADS/THREADS_PER_BLOCK+1;
	size_y_out = ITERATIONS*THREADS*NODES;

	cudaError_t err;
	int device=0;
	cudaDeviceProp device_prop;

	if(err = cudaGetDeviceProperties(&device_prop, device)) cout << cudaGetErrorString(err) << endl;
        cout << "# GPU Device " << device << ": \""<< device_prop.name << "\" with compute capability " << device_prop.major << "." << device_prop.minor << endl;

	double *p, *kij, *Xd, *output;

	if(err = cudaMalloc((void**) &p, NP3*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &kij, NK*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &Xd, size_y_init*sizeof(double))) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMalloc((void**) &output, size_y_out*sizeof(double))) cout << cudaGetErrorString(err) << endl;



	if(err = cudaMemcpy(p, params, NP3*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMemcpy(kij, coupling, NK*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;
	if(err = cudaMemcpy(Xd, y_init, size_y_init*sizeof(double), cudaMemcpyHostToDevice)) cout << cudaGetErrorString(err) << endl;


	cout << "# starting computation with " << BLOCKS << " blocks a " << THREADS_PER_BLOCK << " threads." << endl;
	rk4_three<<<BLOCKS, THREADS_PER_BLOCK>>>(THREADS, Xd, dt, ITERATIONS, stride, p, kij, output);
	if(err = cudaMemcpy(y_out, output, (size_t)size_y_out*sizeof(double), cudaMemcpyDeviceToHost)) cout << "Error in copy device to host: " << cudaGetErrorString(err) << endl;
	

	cudaFree(p);
	cudaFree(kij);
	cudaFree(Xd);
	cudaFree(output);
}



}
