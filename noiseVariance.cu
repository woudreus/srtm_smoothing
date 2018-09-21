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
void NoiseVar(double *dem, double *aggregatedNoiseVar, int *n, int *ncol, int *lengthout);

__global__ void cuImgDiff(const double *dem, double *imgdif, int n, int ncol) {

	int col = blockIdx.x * blockDim.x + (threadIdx.x);
	int row = blockIdx.y * blockDim.y + (threadIdx.y);
	int id  = col + row * ncol;

	int colo = (blockIdx.x * blockDim.x + (threadIdx.x) ) +2;
	int rowo = (blockIdx.y * blockDim.y + (threadIdx.y) ) +2;
	int ido  = colo + rowo * ncol;
	
//	printf("this is cuImgDiff \n\n");

	if ( (col < ncol) && (row < ncol)) {
		double demFocal[] =
		{
			dem[id],
			dem[id + 1],
			dem[id + 2],
			dem[id + 3],
			dem[id + 4],

			dem[id + ncol],
			dem[id + 1 + ncol],
			dem[id + 2 + ncol],
			dem[id + 3 + ncol],
			dem[id + 4 + ncol],

			dem[id + (ncol * 2)],
			dem[id + 1 + (ncol * 2)],
			dem[id + 2 + (ncol * 2)],
			dem[id + 3 + (ncol * 2)],
			dem[id + 4 + (ncol * 2)],

			dem[id + (ncol * 3)],
			dem[id + 1 + (ncol * 3)],
			dem[id + 2 + (ncol * 3)],
			dem[id + 3 + (ncol * 3)],
			dem[id + 4 + (ncol * 3)],

			dem[id + (ncol * 4)],
			dem[id + 1 + (ncol * 4)],
			dem[id + 2 + (ncol * 4)],
			dem[id + 3 + (ncol * 4)],
			dem[id + 4 + (ncol * 4)]
		};

		double focalMat[25] = { 0, 0, 1, 0, 0,
								0, 1, 1, 1, 0,
								1, 1, 0, 1, 1,
								0, 1, 1, 1, 0,
								0, 0, 1, 0, 0 };

		double focalImgSum = 0.0;
		double annulusMat[25];

		for (int i = 0; i < 25; i++) {
			annulusMat[i] = demFocal[i] * focalMat[i];
		}

		for (int i = 0; i < 25; i++) {
			focalImgSum = focalImgSum + annulusMat[i];
		}

//	printf("this is focalImgSum %g \n\n", focalImgSum);
		double imgMean = focalImgSum/25;

	if ((colo > 1 && colo < ncol-1) && (rowo > 1 && rowo < (ncol-1) ))
		imgdif[ido] = (dem[ido] - imgMean);

//	printf("this is imgdif %g \n\n", imgMean);
	}
}

__global__ void cuImgSD(const double *imgdif, double *imgsd, int n, int ncol) {

	int col = blockIdx.x * blockDim.x + (threadIdx.x);
	int row = blockIdx.y * blockDim.y + (threadIdx.y);
	int id = col + row * ncol;

	int colo = (blockIdx.x * blockDim.x + (threadIdx.x ) )+4;
	int rowo = (blockIdx.y * blockDim.y + (threadIdx.y ) )+4;
	int ido = colo + rowo * ncol;

//	printf("this is uImgSD \n\n");
	if ( (col < (ncol-1)) && (col > 1) && (row < (ncol-1)) && (row > 1) ) {
		double difFocal[] =
		{
			imgdif[id],
			imgdif[id + 1],
			imgdif[id + 2],
			imgdif[id + 3],
			imgdif[id + 4],

			imgdif[id + ncol],
			imgdif[id + 1 + ncol],
			imgdif[id + 2 + ncol],
			imgdif[id + 3 + ncol],
			imgdif[id + 4 + ncol],

			imgdif[id + (ncol * 2)],
			imgdif[id + 1 + (ncol * 2)],
			imgdif[id + 2 + (ncol * 2)],
			imgdif[id + 3 + (ncol * 2)],
			imgdif[id + 4 + (ncol * 2)],

			imgdif[id + (ncol * 3)],
			imgdif[id + 1 + (ncol * 3)],
			imgdif[id + 2 + (ncol * 3)],
			imgdif[id + 3 + (ncol * 3)],
			imgdif[id + 4 + (ncol * 3)],

			imgdif[id + (ncol * 4)],
			imgdif[id + 1 + (ncol * 4)],
			imgdif[id + 2 + (ncol * 4)],
			imgdif[id + 3 + (ncol * 4)],
			imgdif[id + 4 + (ncol * 4)]
		};

		double mu = 0.0;
		for (int i = 0; i < 25; i++) 
			mu = mu + difFocal[i]; 
		
		mu = mu / 25;

	//	printf("this is mu %g \n\n",mu);

		double var[25];
		for (int i = 0; i < 25; i++)
			var[i] = pow((difFocal[i] - mu), 2);

		double sumvar = 0.0;
		for (int i = 0; i < 25; i++)
			sumvar = sumvar + var[i];

		double sdv = sqrt(0.04 * sumvar);

		if ((colo > 3 && colo < ncol - 3) && (rowo > 3 && rowo < (ncol - 3)))
		imgsd[ido] = sdv;

//	printf("this is imgsd %g \n\n", sdv);
	}
}
__global__ void cuAgg55(const double *imgSd, double *agg55, int n, int ncol) {
	
	int col = (blockIdx.x * blockDim.x + (threadIdx.x)) * 5;
	int row = (blockIdx.y * blockDim.y + (threadIdx.y)) * 5;
	int id = col + row * ncol;
	
	int colo = (blockIdx.x * blockDim.x + (threadIdx.x));
	int rowo = (blockIdx.y * blockDim.y + (threadIdx.y));
	int outputId = colo + rowo * (ncol/5);

//	printf("this is cuAgg55 \n\n");
	if ( (col > 3) && (col < (ncol - 3) )  && (row > 3) && (row < (ncol - 3) )  ) {	

		double imgSdAgg[] =	
		{
			imgSd[id],
			imgSd[id + 1],
			imgSd[id + 2],
			imgSd[id + 3],
			imgSd[id + 4],

			imgSd[id + ncol],
			imgSd[id + 1 + ncol],
			imgSd[id + 2 + ncol],
			imgSd[id + 3 + ncol],
			imgSd[id + 4 + ncol],

			imgSd[id + (ncol * 2)],
			imgSd[id + 1 + (ncol * 2)],
			imgSd[id + 2 + (ncol * 2)],
			imgSd[id + 3 + (ncol * 2)],
			imgSd[id + 4 + (ncol * 2)],

			imgSd[id + (ncol * 3)],
			imgSd[id + 1 + (ncol * 3)],
			imgSd[id + 2 + (ncol * 3)],
			imgSd[id + 3 + (ncol * 3)],
			imgSd[id + 4 + (ncol * 3)],

			imgSd[id + (ncol * 4)],
			imgSd[id + 1 + (ncol * 4)],
			imgSd[id + 2 + (ncol * 4)],
			imgSd[id + 3 + (ncol * 4)],
			imgSd[id + 4 + (ncol * 4)]
		};

		double sortImgSdAgg[25];
		for (int i = 0; i < 25; i++) {
			int k = 0;
			for (int j = 0; j < 25; j++)
				if (imgSdAgg[i] > imgSdAgg[j])
					k++;
			sortImgSdAgg[k] = imgSdAgg[i];
		}

		double medVal = sortImgSdAgg[13];
	
	if ((colo > 0 && colo < (ncol/5) ) && (rowo > 0 && rowo < (ncol/5) ))
	agg55[outputId] = medVal;
	}
}

__global__ void cuCircFocal(const double *agg55, double *circImgAgg55, int n, int ncol) {

	int col = blockIdx.x * blockDim.x + (threadIdx.x);
	int row = blockIdx.y * blockDim.y + (threadIdx.y);
	int id = col + row * (ncol/5);

	int colo = (blockIdx.x * blockDim.x + (threadIdx.x) )+2;
	int rowo = (blockIdx.y * blockDim.y + (threadIdx.y) )+2;
	int ido = colo + rowo * (ncol/5);

//	printf("this is outputID %i \n\n ", outputId);
//	printf("this is ID %i \n\n ", id);

	if ( (col < (ncol / 5))  && (row < (ncol / 5)) ){
		double focalCirc[] =
		{
			agg55[id],
			agg55[id + 1],
			agg55[id + 2],
			agg55[id + 3],
			agg55[id + 4],

			agg55[id + (ncol / 5)],
			agg55[id + 1 + (ncol / 5)],
			agg55[id + 2 + (ncol / 5)],
			agg55[id + 3 + (ncol / 5)],
			agg55[id + 4 + (ncol / 5)],

			agg55[id + ((ncol / 5) * 2)],
			agg55[id + 1 + ((ncol / 5) * 2)],
			agg55[id + 2 + ((ncol / 5) * 2)],
			agg55[id + 3 + ((ncol / 5) * 2)],
			agg55[id + 4 + ((ncol / 5) * 2)],

			agg55[id + ((ncol / 5) * 3)],
			agg55[id + 1 + ((ncol / 5) * 3)],
			agg55[id + 2 + ((ncol / 5) * 3)],
			agg55[id + 3 + ((ncol / 5) * 3)],
			agg55[id + 4 + ((ncol / 5) * 3)],

			agg55[id + ((ncol / 5) * 4)],
			agg55[id + 1 + ((ncol / 5) * 4)],
			agg55[id + 2 + ((ncol / 5) * 4)],
			agg55[id + 3 + ((ncol / 5) * 4)],
			agg55[id + 4 + ((ncol / 5) * 4)],
		};

		double meanfocal = 0;
		for (int i = 0; i < 25; i++) {
			meanfocal = meanfocal + focalCirc[i];
		}

		meanfocal = meanfocal / 25;


		double medVal = meanfocal;

	// printf("this is medval %g \n\n", medVal);

	if ((colo > 2 && colo < ((ncol/5)-2) ) && (rowo > 2 && rowo < ( (ncol/5) - 2)))
		circImgAgg55[ido] = medVal;

	//printf("this is medval %g \n\n", circImgAgg55[outputId]);
	}
}

using namespace std;
void NoiseVar(double *dem, double *aggregatedNoiseVar, int *n, int *ncol, int *lengthout) {
	
	// initialize device memory variables
	double *d_dem, *d_imgDif, *d_imgSd, *d_agg55, *d_circImgAgg55; //inputs and outputs

	// define grid and total window size
	dim3 grid_in((*ncol / 20), (*ncol / 20)); // grid of 2D blocks
	// these dimensions should be equal to the window dimensions of a raster
	// the product of the 2d window dim should not exceed 1024 threads (depending on your GPU)
	dim3 block_in(20, 20); // block of 2D threads

	// allocate device memory for computations
	// input allocation
	cudaMalloc((void**)&d_dem, *n * sizeof(double));
	// intermediate allocation
	cudaMalloc((void**)&d_imgDif, *n * sizeof(double));
	cudaMalloc((void**)&d_imgSd, *n * sizeof(double));
	cudaMalloc((void**)&d_agg55, *lengthout * sizeof(double));

	// output allocation
	cudaMalloc((void**)&d_circImgAgg55, *lengthout * sizeof(double));

	// copy host memory to allocated device memory
	cudaMemcpy(d_dem, dem, *n * sizeof(double), cudaMemcpyHostToDevice);
	
	// launch kernels with predefined block and thread numbers
	cuImgDiff << <grid_in, block_in >> >(d_dem, d_imgDif, *n, *ncol);
	
	// Free Cuda memory
	cudaFree(d_dem);

	// launch kernel with predefined block and thread numbers
	cuImgSD << <grid_in, block_in >> >(d_imgDif, d_imgSd, *n, *ncol);
	
	// Free Cuda memory
	cudaFree(d_imgDif);

	// launch kernel with predefined block and thread numbers
	cuAgg55 << <grid_in, block_in >> >(d_imgSd, d_agg55, *n, *ncol);

	// Free Cuda memory
	cudaFree(d_imgSd);

	// launch kernel with predefined block and thread numbers
	cuCircFocal << <grid_in, block_in >> >(d_agg55, d_circImgAgg55, *n, *ncol);

//	cudaDeviceSynchronize();

	//perform Array fire opperations
	af::array d_A((*ncol/5),(*ncol/5), d_circImgAgg55, afDevice);

	af::eval(d_A);
	af::sync();

	//interpolation based on specific method
	d_A = resize(5, d_A, AF_INTERP_BILINEAR);

	//return array fire arrays to device memory 
 	double *out = d_A.device<double>();
		
	// copy device memory to host memory
 	cudaMemcpy(aggregatedNoiseVar, out, *n * sizeof(double), cudaMemcpyDeviceToHost);

	// Free Cuda memory
	cudaFree(d_agg55);
	
	cudaFree(d_circImgAgg55);

    cudaFree(out);
}