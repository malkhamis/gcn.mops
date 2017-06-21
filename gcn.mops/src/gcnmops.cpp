// Copyright (C) 2016 Mohammad Alkhamis
// <malkhamis@protonmail.com>

#include <Rmath.h>
#include <Rdefines.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include "macros.h"
#include "gcnmops_priv.h"
#include "gpu/gcnmops_core.h"

#include <unistd.h>

#define SAMPLES_R_DIM 1
#define RANGES_R_DIM 0

// to find how many X_norm ranges can fit in GPU mem
static unsigned int calc_prtn_len(unsigned int nRanges, unsigned int nSamples, unsigned int nClasses, unsigned int gridSize, unsigned int blockSize);

// to package all these args and consts in a single struct
static struct gcnmops_args setup_gcnmops_args (char d_or_h, 
	SEXP x_normS, SEXP IS, SEXP classesS, SEXP covS, SEXP cycS, SEXP idxCN2S, SEXP alphaInitS, SEXP alphaPriorS, SEXP minReadCountS);

using namespace std;


// This is the function called from within R
// x_normS is no longer tranposed
extern "C" SEXP gcnmops_w(SEXP gpuS, SEXP x_normS, SEXP IS, SEXP classesS, SEXP covS, SEXP cycS, SEXP idxCN2S, SEXP alphaInitS,
		SEXP alphaPriorS, SEXP minReadCountS, SEXP returnPosteriorS) {

	unsigned long long i = 0, j = 0; // general purpose

	int dev = *INTEGER(gpuS);
	CHK_CU( cudaSetDevice(dev) );
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, dev);
	Rprintf("Using CUDA-capable Device [%d]: \"%s\"\n", dev, deviceProp.name);
//////////////////////////////////////////////////////////////////////////////////////////////
	// partition the input x_norm according to dev global mem size //
		
	// input GRanges object
	struct X_norm x_norm;
	x_norm.nSamples = INTEGER(GET_DIM(x_normS))[SAMPLES_R_DIM];
	x_norm.nRanges = INTEGER(GET_DIM(x_normS))[RANGES_R_DIM];
	x_norm.data = REAL(x_normS);

	// shortcut since they will be used a lot
	unsigned int N = x_norm.nSamples;
	unsigned int n = length(IS);
	unsigned int nRanges_tot = x_norm.nRanges;
	
	unsigned int x_partn_len = calc_prtn_len(x_norm.nRanges, x_norm.nSamples, n, GRID_SIZE, BLOCK_SIZE);
	struct X_partition_table x_norm_prtnTable =
		_create_partition_table(x_norm.data, x_norm.nRanges, x_norm.nSamples, x_partn_len);
	
	// allocate GPU memory space for an x_norm partition
	struct X_norm d_Xchunk;
	unsigned long long nBytes_xChunk = x_partn_len * x_norm.nSamples * sizeof(double);
	CHK_CU( cudaMalloc((double**)&(d_Xchunk.data), nBytes_xChunk) );
	d_Xchunk.nRanges = x_partn_len;
	d_Xchunk.nSamples = x_norm.nSamples;
	
	Rprintf ("number of kernels to be launched: %d\n", x_norm_prtnTable.nPartitions);
	
//////////////////////////////////////////////////////////////////////////////////////////////
	// unlid and repackage cn.MOPS arguments from SEXP objs

	bool returnPosterior = (bool)(LOGICAL(returnPosteriorS)[0]);
	
	// gcn.mops arguments: allocate mem, package args, and pre-computed consts
	struct gcnmops_args d_gcnmops_args = setup_gcnmops_args ('d', 
		x_normS, IS, classesS, covS, cycS, idxCN2S, alphaInitS, alphaPriorS, minReadCountS);
	struct gcnmops_args h_gcnmops_args = setup_gcnmops_args ('h', 
		x_normS, IS, classesS, covS, cycS, idxCN2S, alphaInitS, alphaPriorS, minReadCountS);

	// copy params to dev's const mem
	cpyArgsToDevSym(h_gcnmops_args);
	
//////////////////////////////////////////////////////////////////////////////////////////////
	// allocate buffers for results and temp vars //
	
	struct gcnmops_buf h_buf = _alloc_gcnmops_buf('h', x_partn_len, n, N);
	struct gcnmops_buf d_buf = _alloc_gcnmops_buf('d', x_partn_len, n, N);

	// dev space for temp vars. Each thread has its own space (assumption: all threads are active and all blocks are resident on SMs)
	struct gcnmops_tmpVars d_workspace = _prealloc_gcnmops_workspace('d', GRID_SIZE*BLOCK_SIZE, n, N);

//////////////////////////////////////////////////////////////////////////////////////////////
	// process the first chunk on device asynchronously
	
	cpyDataToDev(d_Xchunk, x_norm_prtnTable.chunk[0], nRanges_tot);
	printf("processing partition 1 on Device [%d]\n", dev);
	gcnmops_core(GRID_SIZE, BLOCK_SIZE, d_Xchunk, d_gcnmops_args, d_workspace, d_buf);

	// also overlap independent initial host work below with device work
	//////////////////////////////////////////////////////////////////////////////////////////////	
	
	// allocate space for results to be returned to R //
	SEXP RET;
	PROTECT( RET = allocVector(VECSXP, 6) );
	
	SEXP names;
	PROTECT(names = allocVector(STRSXP, 6));
	SET_STRING_ELT(names, 0, mkChar("lambda"));
	SET_STRING_ELT(names, 1, mkChar("alpha"));
	SET_STRING_ELT(names, 2, mkChar("expectedCN"));
	SET_STRING_ELT(names, 3, mkChar("sini"));
	SET_STRING_ELT(names, 4, mkChar("ini"));
	SET_STRING_ELT(names, 5, mkChar("post"));
	setAttrib(RET, R_NamesSymbol, names);
	
	SET_VECTOR_ELT(RET, 0, allocMatrix(REALSXP, x_norm.nRanges, n)); //lambda_est
	SET_VECTOR_ELT(RET, 1, allocMatrix(REALSXP, x_norm.nRanges, n)); //alpha_est
	SET_VECTOR_ELT(RET, 2, allocMatrix(STRSXP , x_norm.nRanges, N)); //expCN
	SET_VECTOR_ELT(RET, 3, allocMatrix(REALSXP, x_norm.nRanges, N)); //ExpLogFoldChange
	SET_VECTOR_ELT(RET, 4, allocVector(REALSXP, x_norm.nRanges   )); //ini


	SEXP postNA;
	PROTECT(postNA = allocVector(LGLSXP, 1));
	LOGICAL(postNA)[0] = R_NaN;
	if(returnPosterior)
		SET_VECTOR_ELT( RET, 5, allocVector(REALSXP, n*N*x_norm.nRanges) ); //alpha_ik = 3D Array (3D-ized at R)
	else
		SET_VECTOR_ELT( RET, 5, postNA ); //alpha_ik = NA

	struct gcnmops_SXP_DATA dataSXP;
	dataSXP.lambda_est = REAL( VECTOR_ELT( RET, 0 ) );
	dataSXP.alpha_est = REAL( VECTOR_ELT( RET, 1 ) );
	dataSXP.ExpLogFoldChange = REAL( VECTOR_ELT( RET, 3 ) );
	dataSXP.ini = REAL( VECTOR_ELT( RET, 4 ) );
	if(returnPosterior) dataSXP.alpha_ik = REAL( VECTOR_ELT( RET, 5 ) );

	//////////////////////////////////////////////////////////////////////////////////////////////
	
//////////////////////////////////////////////////////////////////////////////////////////////
	// start computation iteratively for remaining partitions //
	
	// to convert CN index to CN string  //
	char **CNx = Calloc(n, char*);
	for(i=0; i < n; i++)
		CNx[i] = (char*) CHAR( STRING_ELT(classesS,  i) );
	
	// get a copy so we don't lose them (in case we need them later)
	struct gcnmops_buf h_buf_cpy = h_buf;
	struct gcnmops_SXP_DATA dataSXP_cpy = dataSXP;
	// tracker.. also used to correctly assign CNx strings in RET
	unsigned int nRanges_done = 0;
	
	unsigned int itr;
	for(itr = 0; itr < x_norm_prtnTable.nPartitions; itr++) {
		
		unsigned int prevChunk_len = x_norm_prtnTable.chunk[itr].nRanges;
		
		// if CPU finishes before GPU (barrier)
		devSync();
		
		// result of the previous chunk
		cpyResToHost(h_buf_cpy, d_buf);
		
		// asynchronously, process the next partition (if there is) //
		if((x_norm_prtnTable.nPartitions - itr) > 1) {
			unsigned int curChunk_len = x_norm_prtnTable.chunk[itr+1].nRanges;
			d_Xchunk.nRanges = curChunk_len; // correction for last chunk if it is shorter
			cpyDataToDev(d_Xchunk, x_norm_prtnTable.chunk[itr+1], nRanges_tot);
			printf("processing partition %u on Device [%d]\n", itr+2, dev);
			gcnmops_core(GRID_SIZE, BLOCK_SIZE, d_Xchunk, d_gcnmops_args, d_workspace, d_buf); // non-blocking
		}
		
		// 'if all(x) <= minReadCount' (prev chunk)
		cnmops_leqMRC(x_norm_prtnTable.chunk[itr], h_gcnmops_args, h_buf_cpy, returnPosterior);
		
		// copy results of prev chunk from host buffer to SEXP objects
		//lambda_est
		for(i = 0; i < n; i++) {
			memcpy(dataSXP_cpy.lambda_est, h_buf_cpy.lambda_est, prevChunk_len*sizeof(double));
			dataSXP_cpy.lambda_est += nRanges_tot;
			h_buf_cpy.lambda_est += prevChunk_len;
		}
		//alpha_est
		for(i = 0; i < n; i++) {
			memcpy(dataSXP_cpy.alpha_est, h_buf_cpy.alpha_est, prevChunk_len*sizeof(double));
			dataSXP_cpy.alpha_est += nRanges_tot;
			h_buf_cpy.alpha_est += prevChunk_len;
		}
		//sini
		for(i = 0; i < N; i++) {
			memcpy(dataSXP_cpy.ExpLogFoldChange, h_buf_cpy.ExpLogFoldChange, prevChunk_len*sizeof(double));
			dataSXP_cpy.ExpLogFoldChange += nRanges_tot;
			h_buf_cpy.ExpLogFoldChange += prevChunk_len;
		}
		//alpha_ik
		if(returnPosterior) {
			for(i = 0; i < N*n; i++) {
				memcpy(dataSXP_cpy.alpha_ik, h_buf_cpy.alpha_ik, prevChunk_len*sizeof(double));
				dataSXP_cpy.alpha_ik += nRanges_tot;
				h_buf_cpy.alpha_ik += prevChunk_len;
			}
		}
		// ini
		memcpy(dataSXP_cpy.ini, h_buf_cpy.ini, prevChunk_len*sizeof(double));

		// set expCN string
		unsigned int retCN_idx = nRanges_done;
		for(i=0; i < N; i++) {
			for(j=0; j < prevChunk_len; j++)
				SET_STRING_ELT( VECTOR_ELT(RET,2), retCN_idx + j, mkChar(CNx[ (unsigned int)h_buf_cpy.expCN_idx[i*prevChunk_len+j] ]) );
			retCN_idx += nRanges_tot;
		}

		nRanges_done += prevChunk_len;
		
		// reset pointers
		h_buf_cpy = h_buf;
		// increment SXP data pointers
		dataSXP_cpy.lambda_est = dataSXP.lambda_est + nRanges_done;
		dataSXP_cpy.alpha_est = dataSXP.alpha_est + nRanges_done;
		dataSXP_cpy.ExpLogFoldChange = dataSXP.ExpLogFoldChange + nRanges_done;
		dataSXP_cpy.alpha_ik = dataSXP.alpha_ik + nRanges_done;
		dataSXP_cpy.ini = dataSXP.ini + nRanges_done;
		
	}

//////////////////////////////////////////////////////////////////////////////////////////////
	/* clean up, and exit */
	
	CHK_CU( cudaFree(d_Xchunk.data) );
	Free(CNx);
	_free_gcnmops_buf('h', h_buf);
	_free_gcnmops_buf('d', d_buf);
	_free_gcnmops_args('d', d_gcnmops_args);
	_free_gcnmops_workspace('d', d_workspace);
	
	cudaDeviceReset(); // just in case
	UNPROTECT(3);
	
	return RET;
}


static struct gcnmops_args setup_gcnmops_args (char d_or_h, SEXP x_normS, SEXP IS, SEXP classesS, SEXP covS, SEXP cycS, SEXP idxCN2S, SEXP alphaInitS,
		SEXP alphaPriorS, SEXP minReadCountS) {

	// first, setup args on host //
	
	unsigned int N = INTEGER(GET_DIM(x_normS))[SAMPLES_R_DIM];
	unsigned int n = length(IS);
	double *I=REAL(IS);
	double *alphaInit=REAL(alphaInitS);
	double *alphaPrior=REAL(alphaPriorS);
	double *cov=REAL(covS);
	
	// initially, allocate space on host
	struct gcnmops_args h_args = _alloc_gcnmops_args('h', n, N);
	
	// set constants
	h_args.nRanges_tot = INTEGER(GET_DIM(x_normS))[RANGES_R_DIM];
	h_args.nSamples = N;
	h_args.nClasses = n;
	//h_args.minReadCount = (unsigned int)(REAL(minReadCountS)[0]);
	h_args.minReadCount = REAL(minReadCountS)[0];
	h_args.cyc = INTEGER(cycS)[0];
	h_args.idxCN2 = INTEGER(idxCN2S)[0] - 1;
	
	// cpy/calc arrays
	unsigned int i;
	for(i = 0; i < n; i++) {
		h_args.I[i] = I[i];
		h_args.alphaInit[i] = alphaInit[i];
		h_args.alphaPrior[i] = alphaPrior[i];
		h_args.log2_I[i] = log2(I[i]);	
		h_args.abs_log2_I[i] = fabs(h_args.log2_I[i]);
	}
	for(i = 0; i < N; i++)
		h_args.cov[i] = cov[i];
	
	// then, setup args on device if this request is for device //
	if((d_or_h == 'd') || (d_or_h == 'D')) {

		struct gcnmops_args d_args = _alloc_gcnmops_args('d', n, N);
		d_args.nRanges_tot = h_args.nRanges_tot;
		d_args.nSamples = h_args.nSamples;
		d_args.nClasses = h_args.nClasses;
		d_args.minReadCount = h_args.minReadCount;
		d_args.cyc = h_args.cyc;
		d_args.idxCN2 = h_args.idxCN2;
		
		unsigned int nBytes = N_ELMS_GCNMOPS_ARR_ARGS(N, n) * sizeof(double);
		CHK_CU( cudaMemcpy(d_args.mem_handle, h_args.mem_handle, nBytes, cudaMemcpyHostToDevice) );
	
		_free_gcnmops_args('h', h_args);
		return d_args;
		
	} else
		return h_args;
}

static unsigned int calc_prtn_len(unsigned int nRanges, unsigned int nSamples, unsigned int nClasses, unsigned int gridSize, unsigned int blockSize) {

	// determine the physical memory block size and free memory.
	// physMemBlk_nBytes: to avoid cuda out-of-mem error since 
	// cudaMalloc reserve mem space in blocks of Y bytes.
	const size_t oneMB_nBytes = 1024*1024;
	size_t dev_free_nBytes, physMemBlk_nBytes, temp, unreservable_nBytes; char *dummy;
	CHK_CU( cudaMemGetInfo ( &dev_free_nBytes, NULL ) );
	CHK_CU( cudaMalloc( (char**) &dummy, 1) );
	CHK_CU( cudaMemGetInfo ( &temp, NULL ) );
	cudaFree(dummy);
	//minimum size of physically reservable mem space
	physMemBlk_nBytes = dev_free_nBytes - temp;
	//eliminate unreservable mem space
	unreservable_nBytes = dev_free_nBytes % physMemBlk_nBytes;
	//dev_free_nBytes -= (dev_free_nBytes % physMemBlk_nBytes);
	dev_free_nBytes -= unreservable_nBytes;
#ifdef MAX_DRAM 
	#if MAX_DRAM >= 1
	size_t max_dram_nBytes = ((size_t)MAX_DRAM) * oneMB_nBytes;
	if(dev_free_nBytes > max_dram_nBytes) dev_free_nBytes = max_dram_nBytes;
	#endif
#endif
	if(dev_free_nBytes < oneMB_nBytes) Rprintf("Less then 1MB of device RAM will be used!\n");
	//cudaFree(dummy);
	
	// determine nBytes for each anticipated GPU allocation
	// round nBytes up to make it multiple of 'physMemBlk_nBytes'
    size_t nBytes_X = nSamples * nRanges * sizeof(double);
    size_t nBytesAdj_X = nBytes_X + (physMemBlk_nBytes - (nBytes_X % physMemBlk_nBytes));
    
    size_t nBytes_buf = N_ELMS_GCNMOPS_BUFSLICE(nSamples, nClasses) * nRanges * sizeof(double);
    size_t nBytesAdj_buf = nBytes_buf + (physMemBlk_nBytes - (nBytes_buf % physMemBlk_nBytes));
    
    size_t nBytes_workspace = N_ELMS_GCNMOPS_TMPVARS(nSamples, nClasses) * blockSize * gridSize * sizeof(double);
    nBytes_workspace += (physMemBlk_nBytes - (nBytes_workspace % physMemBlk_nBytes));
    
    size_t nBytes_args = N_ELMS_GCNMOPS_ARR_ARGS(nSamples, nClasses) * sizeof(double);
    nBytes_args += (physMemBlk_nBytes - (nBytes_args % physMemBlk_nBytes));


	size_t fixed_nBytes = nBytes_workspace + nBytes_args;
	size_t var_nBytes = nBytes_X + nBytes_buf;
	size_t var_nBytesAdj = nBytesAdj_X + nBytesAdj_buf;
	size_t nBytes_tot = fixed_nBytes + var_nBytesAdj;

	if(dev_free_nBytes >= nBytes_tot)
		return nRanges;
	else { 
		// 'buf' and 'X' will have separate allocation and potentially reserve extra unused space.
		// Those, we assume that there is free space less 2 times the physMemBlk_nBytes.
		size_t prtn_len = (dev_free_nBytes - (size_t)(fixed_nBytes + 2*physMemBlk_nBytes)) / (size_t)(var_nBytes/nRanges);
#ifdef DEBUG
		if (prtn_len > (size_t)nRanges)
			error("%s:%s @ line %d\nSuggested partition length=%ld was greater than nRanges=%u. This is likely due to insufficient GPU DRAM!\n\n",__FILE__, __func__, __LINE__, prtn_len, nRanges);
#endif
		// errors, if any, will still be caught at '_create_partition_table(..)'
		return( (unsigned int)prtn_len );
	}
}
