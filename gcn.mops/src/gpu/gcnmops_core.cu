/* This file is a substantially-modified version of the 
   following files: cnmops.cpp & cn.mops.R. */

// Copyright (C) 2011 Klambauer Guenter 
// <klambauer@bioinf.jku.at>

// Copyright (C) 2016 Mohammad Alkhamis
// <malkhamis@protonmail.com>

#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include <Rdefines.h>
#include "../gcnmops_typedefs.h"
#include "../macros.h"
#include "math_ud.cuda"
#include <time.h>


__constant__ unsigned int N;
__constant__ unsigned int n;
//__constant__ unsigned int minReadCount;
__constant__ double minReadCount;
__constant__ unsigned int cyc;
__constant__ unsigned int idxCN2;
__constant__ unsigned int nRanges;

__constant__ double magic_val; // 1.0E-10
__constant__ double eps; //1.0E-100

#ifdef USE_CONST_MEM
	#ifndef MAX_nSAMPLES
	#define MAX_nSAMPLES 100
	#endif
	#ifndef MAX_nCLASSES
	#define MAX_nCLASSES 100
	#endif

	__constant__ double I[MAX_nCLASSES];
	__constant__ double alphaInit[MAX_nCLASSES];
	__constant__ double alphaPrior[MAX_nCLASSES];
	__constant__ double log2_I[MAX_nCLASSES];
	__constant__ double abs_log2_I[MAX_nCLASSES];

	__constant__ double cov[MAX_nSAMPLES];
#endif

//#define MAGIC_VAL 1.0E-10


__global__ void _gcnmops_gtMRC(struct X_norm x_chunk, struct gcnmops_args args, struct gcnmops_tmpVars draftspace, struct gcnmops_buf res_buf);

extern "C" void cpyArgsToDevSym(struct gcnmops_args h_args) {

	// constant symbols are defined at the top of this file
	
	CHK_CU( cudaMemcpyToSymbol(N, &h_args.nSamples, sizeof(unsigned int)) );
	CHK_CU( cudaMemcpyToSymbol(n, &h_args.nClasses, sizeof(unsigned int)) );
	//CHK_CU( cudaMemcpyToSymbol(minReadCount, &h_args.minReadCount, sizeof(unsigned int)) );
	CHK_CU( cudaMemcpyToSymbol(minReadCount, &h_args.minReadCount, sizeof(double)) );
	CHK_CU( cudaMemcpyToSymbol(cyc, &h_args.cyc, sizeof(unsigned int)) );
	CHK_CU( cudaMemcpyToSymbol(idxCN2, &h_args.idxCN2, sizeof(unsigned int)) );
	
	double local_eps = 1.0e-100, local_magic_val = 1.0e-10;
	CHK_CU( cudaMemcpyToSymbol(eps, &local_eps, sizeof(double)) );
	CHK_CU( cudaMemcpyToSymbol(magic_val, &local_magic_val, sizeof(double)) );
	
#ifdef USE_CONST_MEM
	if(h_args.nClasses > MAX_nCLASSES) {
		cudaDeviceReset();
		error("number of classes exceeded 'MAX_nCLASSES' as defined. Please, recompile after setting 'MAX_nCALSSES' to a higher value in the Makevar file. Alternatively, you may remove the '-DUSE_CONST_MEM' flag and recompile to let the program use the GPU global memory for storing cn.MOPS arguments and pre-computed constants.\n");
	}
	
	if(h_args.nSamples > MAX_nSAMPLES) {
		cudaDeviceReset();
		error("number of samples exceeded 'MAX_nSAMPLES' as defined. Please, recompile after setting 'MAX_nSAMPLES' to a higher value in the Makevar file. Alternatively, you may remove the '-DUSE_CONST_MEM' flag and recompile to let the program use the GPU global memory for storing cn.MOPS arguments and pre-computed constants.\n");
	}
	
	CHK_CU( cudaMemcpyToSymbol(I, h_args.I, h_args.nClasses*sizeof(double)) );
	CHK_CU( cudaMemcpyToSymbol(alphaInit, h_args.alphaInit, h_args.nClasses*sizeof(double)) );
	CHK_CU( cudaMemcpyToSymbol(alphaPrior, h_args.alphaPrior, h_args.nClasses*sizeof(double)) );
	CHK_CU( cudaMemcpyToSymbol(log2_I, h_args.log2_I, h_args.nClasses*sizeof(double)) );
	CHK_CU( cudaMemcpyToSymbol(abs_log2_I, h_args.abs_log2_I, h_args.nClasses*sizeof(double)) );
	
	CHK_CU( cudaMemcpyToSymbol(cov, h_args.cov, h_args.nSamples*sizeof(double)) );
#endif
	
}


extern "C" void cpyDataToDev(struct X_norm d_Xchunk, struct X_norm h_Xchunk, unsigned int nRanges_tot) {
	// data is stored in a col-major (cols are samples)
	// copy portion of each samples to device memory
	for(unsigned int i = 0; i < h_Xchunk.nSamples; i++) {
		unsigned int h_offset = i * nRanges_tot;
		unsigned int d_offset = i * h_Xchunk.nRanges;
		CHK_CU( cudaMemcpy(d_Xchunk.data + d_offset,   h_Xchunk.data + h_offset,
						   h_Xchunk.nRanges * sizeof(double), cudaMemcpyHostToDevice) );
	}
	
	CHK_CU( cudaMemcpyToSymbol(nRanges, &h_Xchunk.nRanges, sizeof(unsigned int)) );

}


extern "C" void cpyResToHost(struct gcnmops_buf h_buf, struct gcnmops_buf d_buf) {

	CHK_CU( cudaMemcpy(h_buf.mem_handle, d_buf.mem_handle, h_buf.nSlices * h_buf.nBytes_slice, cudaMemcpyDeviceToHost) );
}


extern "C" void devSync() {
	// global barrier
	cudaError_t errAsync = cudaDeviceSynchronize();
	if (errAsync != cudaSuccess) {
		printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));
		cudaDeviceReset();
	}
	//sanity check
	cudaError_t errSync  = cudaGetLastError();
	if (errSync  != cudaSuccess) {
		printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
		cudaDeviceReset();
	}
}

extern "C" void gcnmops_core(unsigned int grid_size, unsigned int block_size,
	struct X_norm d_Xchunk, struct gcnmops_args d_args, struct gcnmops_tmpVars d_workspace, struct gcnmops_buf d_buf) {
		
	CHK_CU( cudaFuncSetCacheConfig(&_gcnmops_gtMRC, cudaFuncCachePreferL1) );	

	// to avoid freezing displays
	//unsigned long long usec = 1000L;
	//struct timespec req = {0};
	//req.tv_sec = 0;
	//req.tv_nsec = 250;
	//nanosleep(&req, (struct timespec *)NULL);
	
	dim3 grid(grid_size);
	dim3 block(block_size);
	_gcnmops_gtMRC <<<grid,block>>> (d_Xchunk, d_args, d_workspace, d_buf);

}

extern "C" void cnmops_leqMRC(struct X_norm h_Xchunk, struct gcnmops_args h_args, struct gcnmops_buf h_buf, bool returnPosterior) {

	bool x_leq_mrc;
	for(unsigned int row = 0; row < h_Xchunk.nRanges; row++) {
	
		x_leq_mrc = true;
		for(unsigned int col = 0; col < h_Xchunk.nSamples; col++) {
			if(h_Xchunk.data[row + col*h_args.nRanges_tot] > h_args.minReadCount) {
				x_leq_mrc = false;
				break;
			}
		}
		
		if(x_leq_mrc == true) {
			set_arrVal(h_buf.lambda_est+row, h_args.nClasses, h_Xchunk.nRanges, 0);
	
			set_arrVal(h_buf.alpha_est+row, h_args.nClasses, h_Xchunk.nRanges, 0);
			h_buf.alpha_est[row + h_args.idxCN2*h_Xchunk.nRanges] = 1;

			// assign the CN index. Make strings at the wrapper
			set_arrVal(h_buf.expCN_idx+row, h_args.nSamples, h_Xchunk.nRanges, h_args.idxCN2);
			
			h_buf.ini[row] = 0.0;
	
			set_arrVal(h_buf.ExpLogFoldChange+row, h_args.nSamples, h_Xchunk.nRanges, 0);
	
			// set only idxCN2'th row of posterior to 1
			if(returnPosterior) {
				set_arrVal(h_buf.alpha_ik + row, h_args.nClasses*h_args.nSamples, h_Xchunk.nRanges, 0);
				set_arrVal(h_buf.alpha_ik + row + h_args.idxCN2*h_Xchunk.nRanges, h_args.nSamples, h_args.nClasses*h_Xchunk.nRanges, 1);
			}
		}
	}
	
}


__global__ void _gcnmops_gtMRC(struct X_norm x_chunk, struct gcnmops_args args, struct gcnmops_tmpVars draftspace, struct gcnmops_buf res_buf) {

	unsigned int g_tid = blockIdx.x * blockDim.x + threadIdx.x;
	if(!(g_tid < nRanges)) return;
	
#ifndef USE_CONST_MEM
	double *I = args.I;
	double *alphaInit = args.alphaInit;
	double *alphaPrior = args.alphaPrior;
	double *log2_I = args.log2_I;
	double *abs_log2_I = args.abs_log2_I;
	double *cov = args.cov;
#endif
	
	struct gcnmops_buf local_buf = res_buf;
	
	double *x = x_chunk.data;

	double *x_over_cov_sorted = draftspace.x_over_cov_sorted + g_tid;
	double *lambdaInit = draftspace.lambdaInit + g_tid;
	double *lg = draftspace.lg + g_tid;
	double *alpha_i = draftspace.alpha_i + g_tid;
	double *draftspace_ini = draftspace.draftspace_ini + g_tid;
	

	// intermediates
	//double eps=1e-100;
	int i, k, idx;
	
	// map input and output to threads
	x += g_tid;
	INC_GCNMOPS_D_BUF(local_buf, g_tid);
	
	// grid-strided loop in case there is more rows than total threads
	unsigned int nThrsPerGrid = gridDim.x * blockDim.x;
	unsigned long long row;
	for(row = g_tid; row < nRanges; row += nThrsPerGrid) {

			double sumx=0.0;
	
			for(k = 0; k < N; k++) {
				lg[k*nThrsPerGrid] = lgamma(x[k*nRanges]+1);
				sumx += x[k*nRanges];
			}
			
			double sumAlphaPrior=0.0;
			for(i = 0; i < n; i++)
				sumAlphaPrior=sumAlphaPrior+alphaPrior[i];

			for(i = 0; i < n; i++)
				local_buf.alpha_est[i*nRanges]=alphaInit[i];

			
			// calculate lambdaInit /////////////////////////////////////////////////////
			for(i = 0; i < N; i++)
				 x_over_cov_sorted[i*nThrsPerGrid] = x[i*nRanges] / cov[i];
			// sort the data to find the median
			bubsort_inPlace (x_over_cov_sorted, N, nThrsPerGrid);
			double lambda = median_ud(x_over_cov_sorted, N, nThrsPerGrid);
			if (lambda < magic_val)
				lambda = max(mean_ud(x_over_cov_sorted, N, nThrsPerGrid), 1.0);
			mul_const_vec(I, n, lambda, lambdaInit, nThrsPerGrid);
			/////////////////////////////////////////////////////////////////////////////
			
			for(i = 0; i < n; i++)
				local_buf.lambda_est[i*nRanges]=lambdaInit[i*nThrsPerGrid]*I[i];
			
			double colSum=0.0;
			double sumIalpha_i=0.0;
			for(int cycNr=0; cycNr < cyc; cycNr++) {
				for (k=0; k<N; k++) { //next col
					colSum=0;
					for (i=0; i<n; i++){
						//next elm in current col
						local_buf.alpha_ik[(k*n+i)*nRanges] = local_buf.alpha_est[i*nRanges]*exp(
																				(x[k*nRanges]*log(cov[k]*local_buf.lambda_est[i*nRanges])) 
																				- lg[k*nThrsPerGrid] 
																				- cov[k]*local_buf.lambda_est[i*nRanges] );

						if (local_buf.alpha_ik[(k*n+i)*nRanges] < eps) local_buf.alpha_ik[(k*n+i)*nRanges] = eps;

						colSum=colSum+local_buf.alpha_ik[(k*n+i)*nRanges];
					}
					for (i=0; i<n; i++)
						local_buf.alpha_ik[(k*n+i)*nRanges]=local_buf.alpha_ik[(k*n+i)*nRanges]/(colSum);
				}

				sumIalpha_i=0.0;
				for (i=0; i<n; i++) {
					alpha_i[i*nThrsPerGrid]=0.0;
					for (k=0; k<N; k++){
						sumIalpha_i=sumIalpha_i+I[i]*local_buf.alpha_ik[(k*n+i)*nRanges]*cov[k];
						alpha_i[i*nThrsPerGrid]=alpha_i[i*nThrsPerGrid]+local_buf.alpha_ik[(k*n+i)*nRanges];
					}
					alpha_i[i*nThrsPerGrid]=alpha_i[i*nThrsPerGrid]/N;
				}

				for(i=0; i<n; i++)
					local_buf.alpha_est[i*nRanges]=(alphaPrior[i]+alpha_i[i*nThrsPerGrid])/(1+sumAlphaPrior);
				for(i=0; i<n; i++)
					local_buf.lambda_est[i*nRanges]=sumx/sumIalpha_i*I[i];
			}
		
			// calculate ini ////////////////////////////////////////////////////
			mat_mul_R(abs_log2_I, 1, n, 1, local_buf.alpha_ik, n, N, nRanges, draftspace_ini, nThrsPerGrid);
			*(local_buf.ini) = mean_ud(draftspace_ini, N, nThrsPerGrid);	
			/////////////////////////////////////////////////////////////////////
			
			// calculate ExpLogFoldChange ///////////////////////////////////////
			mat_mul_R(log2_I, 1, n, 1, local_buf.alpha_ik, n, N, nRanges, local_buf.ExpLogFoldChange, nRanges);
			/////////////////////////////////////////////////////////////////////
			
			// assign the CN index. Make strings at the wrapper /////////////////
			for(idx=0; idx < N; idx++)
				local_buf.expCN_idx[idx*nRanges] = get_idxOfMax_col(local_buf.alpha_ik+(idx*n*nRanges), n, nRanges);
			/////////////////////////////////////////////////////////////////////
//		}
	
		/* increment buffer and point to the next output for this chunk */
		INC_GCNMOPS_D_BUF(local_buf, nThrsPerGrid);
		x += nThrsPerGrid;
	}

	//Free(alpha_i); // in original cn.mops, this is returned, but never used!!
}

