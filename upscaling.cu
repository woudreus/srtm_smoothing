#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <device_launch_parameters.h>
#include <arrayfire.h>
#include <af/cuda.h>

extern "C"
__declspec(dllexport)
void bicubUpscaling(double *x,	double *outx, int *ncol, int *n);

void bicubUpscaling(double *x, double *outx, int *ncol, int *n) {
	// initialize device memory variables
	double *d_x;

	//input memory allocation
	cudaMalloc((void**)&d_x, (*n / 9) * sizeof(double));

	//copy host memory to allocated device memory
	cudaMemcpy(d_x, x, (*n / 9) * sizeof(double), cudaMemcpyHostToDevice);

	//conduct AF operations
	af::array d_A((*ncol / 3), (*ncol / 3), d_x, afDevice);
	af::eval(d_A);
	af::sync();
	d_A = resize(3, d_A, AF_INTERP_BICUBIC);

	//return AF arrays to device memory 
	double *x_interp = d_A.device<double>();

	//copy device memory to host memory
	cudaMemcpy(outx, x_interp, *n * sizeof(double), cudaMemcpyDeviceToHost);

	cudaFree(d_x);
	cudaFree(x_interp);

}
