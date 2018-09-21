#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <device_launch_parameters.h>
#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>
#include <sstream>


#define CUDA_ERROR_CHECK
#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

extern "C" 
__declspec(dllexport)
void RCUDA2(double *z, double *noise, double *noisesq, double *variance, double *n, int *count, int *ncol, int *n_out,
			double *outz, double *outw, double *outwsq, double *outn, double *outvg, double *outneff, double *outvm,
			double *outvstat);


__global__ void gallantSmoothing(const double *A, const double *B, const double *C, const double *D, const double *E,
	double *outz, double *outw, double *outwsq, double *outvg, double *outn,
	int x, int ncol, double *outvm, double *outneff, double *outvstat) {


	int col = (blockIdx.x * blockDim.x + threadIdx.x) * 3;
	int row = (blockIdx.y * blockDim.y + threadIdx.y) * 3;
	int id = col + row * (ncol);

	int colo = (blockIdx.x * blockDim.x + (threadIdx.x));
	int rowo = (blockIdx.y * blockDim.y + (threadIdx.y));
	int ido = colo + rowo * (ncol / 3);

	//input id

	 //	printf("this is the id %d \n", id);

	// define focal window processing dimensions

	if ( ( (colo > 0) && (colo < (ncol / 3) )) && ( (rowo > 0) && (rowo < (ncol / 3) ))) {

		double zbarWin[] =
		{
			A[id],
			A[id + 1],
			A[id + 2],

			A[id + ncol],
			A[id + 1 + ncol],
			A[id + 2 + ncol],

			A[id + (ncol * 2)],
			A[id + 1 + (ncol * 2)],
			A[id + 2 + (ncol * 2)]
		};

		double wWin[] =
		{
			B[id],
			B[id + 1],
			B[id + 2],

			B[id + ncol],
			B[id + 1 + ncol],
			B[id + 2 + ncol],

			B[id + (ncol * 2)],
			B[id + 1 + (ncol * 2)],
			B[id + 2 + (ncol * 2)]
		};

		double wsqWin[] =
		{
			C[id],
			C[id + 1],
			C[id + 2],

			C[id + ncol],
			C[id + 1 + ncol],
			C[id + 2 + ncol],

			C[id + (ncol * 2)],
			C[id + 1 + (ncol * 2)],
			C[id + 2 + (ncol * 2)]
		};

		double vgWin[] =
		{
			D[id],
			D[id + 1],
			D[id + 2],

			D[id + ncol],
			D[id + 1 + ncol],
			D[id + 2 + ncol],

			D[id + (ncol * 2)],
			D[id + 1 + (ncol * 2)],
			D[id + 2 + (ncol * 2)]
		};

		double nWin[] =
		{
			E[id],
			E[id + 1],
			E[id + 2],

			E[id + ncol],
			E[id + 1 + ncol],
			E[id + 2 + ncol],

			E[id + (ncol * 2)],
			E[id + 1 + (ncol * 2)],
			E[id + 2 + (ncol * 2)]
		};

		//  printf("this is wwin: %g \n", double(wWin[0]) );

		// execute w and wsq
		double w = 0.0;
		for (int i = 0; i < 9; i++) {
			w = w + wWin[i];
		}


		outw[ido] = w;


	//		printf("this is w: %g \n", w);

		double wsq = 0.0;
		for (int i = 0; i < 9; i++) {
			wsq = wsq + wsqWin[i];
		}

			//	if ((colo > 0 && colo < (ncol / 3)) && (rowo > 0 && rowo < (ncol / 3))) {
			outwsq[ido] = wsq;
			//	}

				//	printf("this is wsq: %d \n", wsq);
				//  execute zbar

			double tempZbar = 0.0;
			double tempWZ[9];

			for (int i = 0; i < 9; i++) {
				tempWZ[i] = (zbarWin[i] * wWin[i])/w;
			}

			for (int i = 0; i < 9; i++) {
				tempZbar = tempZbar + tempWZ[i];
			}

			double zbar = tempZbar ;

			//	printf("this is zbar: %g \n", zbar);

		//	if ((colo > 0 && colo < (ncol / 3)) && (rowo > 0 && rowo < (ncol / 3))) {
			outz[ido] = zbar;
			//	}

				//  execute vbg
			double tempvbgArray[9];
			for (int i = 0; i < 9; i++) {

				double zbarDif= (zbarWin[i] - zbar);

				tempvbgArray[i] = (wWin[i] * pow(zbarDif, 2) )/w;
			}
			double tempvbgNum = 0.0;
			for (int i = 0; i < 9; i++) {
				tempvbgNum = (tempvbgNum + tempvbgArray[i]);
			}
			double vbg = (tempvbgNum );

		/*	if (vbg > 20) {
				printf("this is vbg: %g \n", vbg);
			}
			*/
			// execute vwg
			double tempvwgArray[9];
			for (int i = 0; i < 9; i++) {
				tempvwgArray[i] = (wWin[i] * vgWin[i]);
			}

			double tempvwg = 0;

			for (int i = 0; i < 9; i++) {
				tempvwg = (tempvwg + tempvwgArray[i])/w;
			}

			double vwg = (tempvwg);

			//	printf("this is vwg: %g \n", vwg);

			// execute vg

			double vg = (vbg + vwg);

			//	if ((colo > 0 && colo < (ncol / 3)) && (rowo > 0 && rowo < (ncol / 3))) {
			outvg[ido] = vg;
			//	}
				//		printf("this is vg: %g \n", vg);

				// execute vm
			//	if ((colo > 0 && colo < (ncol / 3)) && (rowo > 0 && rowo < (ncol / 3))) {
			outvm[ido] = 1 / w;
			//	}

				//		printf("this is vm: %d \n", vm);

				// execute n

			double n = 0.0;
			for (int i = 0; i < 9; i++) {
				n = n + nWin[i];
			}

			//	if ((colo > 0 && colo < (ncol / 3)) && (rowo > 0 && rowo < (ncol / 3))) {
			outn[ido] = n;
			//	}

				// execute neff
			//	if ((colo > 0 && colo < (ncol / 3)) && (rowo > 0 && rowo < (ncol / 3))) {
			outneff[ido] = (pow(w, 2) / wsq);
			//	}

					//		printf("this is neff: %g \n", neff);

					// execute mv

			double mv = n / w;
			//		printf("this is mv: %g \n", mv);

			// get vstats
	//		if ((colo > 0 && colo < (ncol / 3)) && (rowo > 0 && rowo < (ncol / 3))) {
			outvstat[ido] = vg / mv;
		}
	}
	

	inline void __cudaSafeCall(cudaError err, const char *file, const int line)
	{
#ifdef CUDA_ERROR_CHECK
		if (cudaSuccess != err)
		{
			fprintf(stderr, "cudaSafeCall() failed at %s:%i : %s\n",
				file, line, cudaGetErrorString(err));
		//	exit(-1);
		}
#endif

		return;
	}

	inline void __cudaCheckError(const char *file, const int line)
	{
#ifdef CUDA_ERROR_CHECK
		cudaError err = cudaGetLastError();
		if (cudaSuccess != err)
		{
			fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
				file, line, cudaGetErrorString(err));
		//	exit(-1);
		}

		// More careful checking. However, this will affect performance.
		// Comment away if needed.
		err = cudaDeviceSynchronize();
		if (cudaSuccess != err)
		{
			fprintf(stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
				file, line, cudaGetErrorString(err));
//			exit(-1);
		}

#endif


		return;
	}

 // call host function 
void RCUDA2(double *z, double *noise, double *noisesq, double *variance, double *n, int *count, int *ncol, int *n_out,
	double *outz, double *outw, double *outwsq, double *outn, double *outvg, double *outneff, double *outvm,
	double *outvstat) {


	// initialize device memory variables
	double *d_n, *d_vg, *d_wsq, *d_z, *d_w, *d_outz, *d_outw, *d_outwsq, *d_outvg, *d_outn,
		   *d_vm, *d_neff, *d_vstat, *d_xcrit; //inputs and outputs

	//define grid and total window size
	dim3 grid((*ncol / 24), (*ncol / 24)); // grid of 2D blocks

	//these dimensions should be equal to the window dimensions of a raster
	// the product of the 2d window dim should not exceed 1024 threads (depending on your GPU)
	dim3 block(24, 24); // block of 2D threads

	// allocate device memory for computations
	// input allocation
	CudaSafeCall(cudaMalloc((void**)&d_z, *count * sizeof(double)));
	
	CudaSafeCall(cudaMalloc((void**)&d_w, *count * sizeof(double)));
	CudaSafeCall(cudaMalloc((void**)&d_wsq, *count * sizeof(double)));
	CudaSafeCall(cudaMalloc((void**)&d_vg, *count * sizeof(double)));
	CudaSafeCall(cudaMalloc((void**)&d_n, *count * sizeof(double)));
	CudaSafeCall(cudaMalloc((void**)&d_xcrit, *count * sizeof(double)));

	//intermediate allocations
	CudaSafeCall(cudaMalloc((void**)&d_vm, *n_out * sizeof(double)));
	CudaSafeCall(cudaMalloc((void**)&d_neff, *n_out * sizeof(double)));
	CudaSafeCall(cudaMalloc((void**)&d_vstat, *n_out * sizeof(double)));

	// output allocation
	CudaSafeCall(cudaMalloc((void**)&d_outz, *n_out * sizeof(double)));
	CudaSafeCall(cudaMalloc((void**)&d_outw, *n_out * sizeof(double)));
	CudaSafeCall(cudaMalloc((void**)&d_outwsq, *n_out * sizeof(double)));
	CudaSafeCall(cudaMalloc((void**)&d_outvg, *n_out * sizeof(double)));
	CudaSafeCall(cudaMalloc((void**)&d_outn, *n_out * sizeof(double)));
	
	// copy host memory to allocated device memory
	CudaSafeCall(cudaMemcpy(d_z, z, *count * sizeof(double), cudaMemcpyHostToDevice));
	CudaSafeCall(cudaMemcpy(d_w, noise, *count * sizeof(double), cudaMemcpyHostToDevice));
	CudaSafeCall(cudaMemcpy(d_wsq, noisesq, *count * sizeof(double), cudaMemcpyHostToDevice));
	CudaSafeCall(cudaMemcpy(d_vg, variance, *count * sizeof(double), cudaMemcpyHostToDevice));
	CudaSafeCall(cudaMemcpy(d_n, n, *count * sizeof(double), cudaMemcpyHostToDevice));
	


	
//	printf("this string works");
	
		// call some_function which may throw something
		// launch kernel with predefined block and thread numbers
	gallantSmoothing << <grid, block >> > (d_z, d_w, d_wsq, d_vg, d_n, d_outz, d_outw, d_outwsq, d_outvg, d_outn, *count, *ncol,
			d_vm, d_neff, d_vstat);
	cudaDeviceSynchronize();

//	gallantSmoothing << <grid, block >> > (d_z, *ncol);
		CudaCheckError();
	
	
	// copy device memory to host memory
		CudaSafeCall(cudaMemcpy(outz, d_outz, *n_out * sizeof(double), cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(outw, d_outw, *n_out * sizeof(double), cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(outwsq, d_outwsq, *n_out * sizeof(double), cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(outvg, d_outvg, *n_out * sizeof(double), cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(outn, d_outn, *n_out * sizeof(double), cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(outneff, d_neff, *n_out * sizeof(double), cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(outvm, d_vm, *n_out * sizeof(double), cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(outvstat, d_vstat, *n_out * sizeof(double), cudaMemcpyDeviceToHost));

	cudaFree(d_z);
	cudaFree(d_w);
	cudaFree(d_wsq);
	cudaFree(d_vg);
	cudaFree(d_n);
	cudaFree(d_outz);
	cudaFree(d_outw);
	cudaFree(d_outwsq);
	cudaFree(d_xcrit);
	cudaFree(d_vm);
	cudaFree(d_vstat);
	cudaFree(d_outvg);
	cudaFree(d_outn);
	cudaFree(d_neff);
	
	cudaThreadExit();
}