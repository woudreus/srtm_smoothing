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
void treeOffset(double *dtm, double *palshh, double *hans, double *nonf, double *trh, double *dem,
	int *ncol, int *n);
	//, double *phv_test, double *hans_test, double *blur_test, double *combdata_test);

using namespace af;

__global__ void cuLinearTransformation(const double *cupalshh, double *cuphh, int ncol) {

//	printf("This works \n\n");
	int col = blockIdx.x * blockDim.x + (threadIdx.x);
	int row = blockIdx.y * blockDim.y + (threadIdx.y);
	int id = col + row * ncol;
	
	if (col < ncol-1 && row < ncol-1) {

		if (cupalshh[id] < 6.9) {

			cuphh[id] = 0;

		}
		else if (cupalshh[id] >= 6.9 && cupalshh[id] <= 7.6) {

			cuphh[id] = (cupalshh[id] - 6.9) / 0.7;

		}
		else if (cupalshh[id] > 7.6) {

			cuphh[id] = 1;

		}
	}
}

__global__ void cuHansTransformation(const double *cuhans, double *culinhans, int ncol) {

	int col = blockIdx.x * blockDim.x + (threadIdx.x);
	int row = blockIdx.y * blockDim.y + (threadIdx.y);
	int id = col + row * ncol;

	if (col < ncol-1 && row < ncol-1) {

		if (cuhans[id] <= 50) {

			culinhans[id] = 0;
		}else if (cuhans[id] > 50 && cuhans[id] < 70) {

			culinhans[id] = (cuhans[id] - 50) / 20;
		}else if (cuhans[id] >= 70) {

			culinhans[id] = 1;

		}
	}
}

__global__ void cuLinearFragmentation(const double *cupalshh, double *culinfragHV, int ncol) {

	int col = blockIdx.x * blockDim.x + (threadIdx.x);
	int row = blockIdx.y * blockDim.y + (threadIdx.y);
	int id = col + row * ncol;

	int colo = (blockIdx.x * blockDim.x + (threadIdx.x)) + 2;
	int rowo = (blockIdx.y * blockDim.y + (threadIdx.y)) + 2;
	int ido = colo + rowo * ncol;

	if (col < ncol-1 && row < ncol-1) {
		double focalpalshh[] =
		{
			cupalshh[id],
			cupalshh[id + 1],
			cupalshh[id + 2],
			cupalshh[id + 3],
			cupalshh[id + 4],

			cupalshh[id + ncol],
			cupalshh[id + 1 + ncol],
			cupalshh[id + 2 + ncol],
			cupalshh[id + 3 + ncol],
			cupalshh[id + 4 + ncol],

			cupalshh[id + (ncol * 2)],
			cupalshh[id + 1 + (ncol * 2)],
			cupalshh[id + 2 + (ncol * 2)],
			cupalshh[id + 3 + (ncol * 2)],
			cupalshh[id + 4 + (ncol * 2)],

			cupalshh[id + (ncol * 3)],
			cupalshh[id + 1 + (ncol * 3)],
			cupalshh[id + 2 + (ncol * 3)],
			cupalshh[id + 3 + (ncol * 3)],
			cupalshh[id + 4 + (ncol * 3)],

			cupalshh[id + (ncol * 4)],
			cupalshh[id + 1 + (ncol * 4)],
			cupalshh[id + 2 + (ncol * 4)],
			cupalshh[id + 3 + (ncol * 4)],
			cupalshh[id + 4 + (ncol * 4)]
		};

		double focalV[25];
		double focalH[25];

		double cuverFun[25] ={
		-1.937, 0.801, 2.274, 0.801, -1.937,
		-1.937, 0.801, 2.274, 0.801, -1.937,
		-1.937, 0.801, 2.274, 0.801, -1.937,
		-1.937, 0.801, 2.274, 0.801, -1.937,
		-1.937, 0.801, 2.274, 0.801, -1.937 };

		double cuhorFun[25] = {
		-1.937, -1.937, -1.937, -1.937, -1.937,
		 0.801,  0.801,  0.801,  0.801,  0.801,
		 2.274,  2.274,  2.274,  2.274,  2.274,
		 0.801,  0.801,  0.801,  0.801,  0.801,
	    -1.937, -1.937, -1.937, -1.937, -1.937 };
		

	for (int i = 0; i < 25; i++) {
		focalV[i] = focalpalshh[i] * cuverFun[i];
		focalH[i] = focalpalshh[i] * cuhorFun[i];
	}

	double sumV;
	double sumH;

	for (int i = 0; i < 25; i++) {
		sumV = focalV[i] + sumV;
		sumH = focalH[i] + sumH;
	}

	double linfragHV = sumH + sumV;

	if ( (colo > 1 && colo < ncol - 1 ) && ( rowo > 1 && rowo < (ncol-1) ) ){

	if (linfragHV < -1.5) {

		culinfragHV[ido] = 0;
	}
	else if (linfragHV >= -1.5 && linfragHV <= 0) {

		culinfragHV[ido] = linfragHV/-1.5;

	}
	else if (linfragHV > 0) {

		culinfragHV[ido] =1;
	
	}
  }
}
}
	__global__ void cuCombData(const double *cuphh, const double *culinfraghv, const double *culinhans, const double *cunonf, double *cucombdat, int ncol) {

		int col = blockIdx.x * blockDim.x + (threadIdx.x);
		int row = blockIdx.y * blockDim.y + (threadIdx.y);
		int id = col + row * ncol;

		if (col < ncol-1 && row < ncol-1) {

		double sumCombDat = cuphh[id] + culinfraghv[id] + culinhans[id] + (cunonf[id]);

		if (sumCombDat <=2.65) {

			cucombdat[id] = 0;

		}else if (sumCombDat > 2.5) {

			cucombdat[id] = 1;
		}

	//	printf("this is the cucombdat %g  \n\n", cucombdat[id]);

		}
	}

	af::array gaussianblur(const af::array &in, int window_width, int window_height, int sigma) {
		af::array g = af::gaussiankernel(window_width, window_height, sigma, sigma);
		return convolve(in, g);
		}


	__global__ void cuDEM(const double *cucombdata, const double *cutrh, const double *cudtm, double *dem, int ncol) {
	
		int col = blockIdx.x * blockDim.x + threadIdx.x;
		int row = blockIdx.y * blockDim.y + threadIdx.y;
		int id = col + row * ncol;

		if ((col < (ncol - 1) && (col > 1)) && (row < (ncol - 1) && (row > 1)) ) {
			double treeHeightModel = cucombdata[id] * cutrh[id];
			dem[id] = cudtm[id] - (0.5* treeHeightModel);
	}
}

void treeOffset(double *dtm, double *palshh, double *hans, double *nonf, double *trh, double *dem,
	int *ncol, int *n){
//, double *phv_test, double *hans_test, double *blur_test, double *combdata_test ) {

	// initialize device memory variables
	double *d_dtm, *d_palshh, *d_hans, *d_nonf, *d_trh, *d_dem, *d_linfragHV, *d_phh, *d_linhans, *d_combdata;  //inputs and outputs

																															 //define grid and total window size
	dim3 grid(*ncol/20, *ncol/20); // grid of 2D blocks

	//these dimensions should be equal to the window dimensions of a raster
	// the product of the 2d window dim should not exceed 1024 threads (depending on your GPU)
	dim3 block(20, 20);// block of 2D threads

	// allocate device memory for computations
	// input allocation
	cudaMalloc((void**)&d_dtm, *n * sizeof(double));
	cudaMalloc((void**)&d_palshh, *n * sizeof(double));
	cudaMalloc((void**)&d_hans, *n * sizeof(double));
	cudaMalloc((void**)&d_nonf, *n * sizeof(double));
	cudaMalloc((void**)&d_trh, *n * sizeof(double));

	// intermediary allocations
	cudaMalloc((void**)&d_linfragHV, *n * sizeof(double));
	cudaMalloc((void**)&d_linhans, *n * sizeof(double));
	cudaMalloc((void**)&d_combdata, *n * sizeof(double));
	cudaMalloc((void**)&d_phh, *n * sizeof(double));
	// output allocation
	cudaMalloc((void**)&d_dem, *n * sizeof(double));

	// copy host memory to allocated device memory
	cudaMemcpy(d_dtm, dtm, *n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_palshh, palshh, *n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_hans, hans, *n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_nonf, nonf, *n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_trh, trh, *n * sizeof(double), cudaMemcpyHostToDevice);
	
	// launch kernel with predefined block and thread numbers
	cuLinearTransformation << <grid, block >> >(d_palshh, d_phh, *ncol);

	// launch kernel with predefined block and thread numbers
	cuHansTransformation << <grid, block >> >(d_hans, d_linhans, *ncol);

	// launch kernel with predefined block and thread numbers
	cuLinearFragmentation << <grid, block >> >(d_palshh, d_linfragHV, *ncol);

	// launch kernel with predefined block and thread numbers
	cuCombData << <grid, block >> >(d_phh, d_linfragHV, d_linhans, d_nonf, d_combdata, *ncol);

	cudaFree(d_nonf);
	cudaFree(d_linhans);
	cudaFree(d_linfragHV);
	cudaFree(d_palshh);
	cudaFree(d_phh);
	cudaFree(d_hans);

	//conduct array fire operations
	af::array d_A(*ncol, *ncol, d_combdata, afDevice);

	af::eval(d_A);
	af::sync();

	af::array d_B = gaussianblur(d_A, 5, 5, 1.5);

	//return array fire arrays to device memory 
	double *d_blurfunction = d_B.device<double>();

	// launch kernel with predefined block and thread numbers
	cuDEM << <grid, block >> >(d_blurfunction, d_trh, d_dtm, d_dem, *ncol);

	cudaMemcpy(dem, d_dem, *n * sizeof(double), cudaMemcpyDeviceToHost);
/*	cudaMemcpy(phv_test, d_linfragHV, *n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(hans_test, d_linhans, *n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(combdata_test, d_combdata, *n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(blur_test, d_blurdata, *n * sizeof(double), cudaMemcpyDeviceToHost);
*/
	// free memory
	cudaFree(d_dem);
	cudaFree(d_dtm);
	cudaFree(d_trh);
	cudaFree(d_combdata);
	cudaFree(d_blurfunction);

}