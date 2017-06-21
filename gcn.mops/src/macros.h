#include <sys/time.h>
#include <Rdefines.h>

#ifndef _UD_MACROS_H
#define _UD_MACROS_H

#define XSTR(x) STR(x)
#define STR(x) #x

// specific for nVidia Tesla C2070 (temporary solution)
#ifndef BLOCK_SIZE
	#define BLOCK_SIZE 64
#endif
#ifndef GRID_SIZE
	#define NUM_SM 14
	#define MAX_THREADS_PER_SM 1536
	#define GRID_SIZE ((MAX_THREADS_PER_SM / BLOCK_SIZE) * NUM_SM) // (1536 / 64) * 14 = 336
#endif

#define CHK_CU(call) \
{ \
    const cudaError_t cu_error = call; \
    if (cu_error != cudaSuccess) \
    { \
    	cudaDeviceReset(); \
        error("%s:%s @ line %d\ncode: %d, reason: %s\n", __FILE__, __func__, __LINE__, cu_error, cudaGetErrorString(cu_error));\
    } \
}

// len of expCN, ExpLog = nSamples; len of lambda_est, alpha_est = nClasses, len of: alpha_ik = nSamples*nClasses, len of ini = 1
#define N_ELMS_GCNMOPS_BUFSLICE(nSamples, nClasses) (2*nClasses + 2*nSamples + 1 + nClasses*nSamples)
// len of: cov, I, log2_I, abs_log2_I, alphaPrior = nClasses ; len of alphaInit = nSamples
#define N_ELMS_GCNMOPS_ARR_ARGS(nSamples, nClasses) (nSamples + 5*nClasses)
// len of: lambdaInit (nClasses), x_over_cov_sortede (nSamples), lg (nSamples), alpha_i (nClasses), draftspace_ini (nSamples)
#define N_ELMS_GCNMOPS_TMPVARS(nSamples, nClasses) (3*nSamples + 2*nClasses)


#define INC_GCNMOPS_H_BUF(buffer, inc_val, numClasses, numSamples) \
{ \
	unsigned int MACROVAR__n_inc_val = inc_val * numClasses; \
	unsigned int MACROVAR__N_inc_val = inc_val * numSamples; \
	unsigned int MACROVAR__nN_inc_val = inc_val * numClasses * numSamples; \
\
	buffer.lambda_est += MACROVAR__n_inc_val; \
	buffer.alpha_est += MACROVAR__n_inc_val; \
	buffer.expCN_idx += MACROVAR__N_inc_val; \
	buffer.ExpLogFoldChange += MACROVAR__N_inc_val; \
	buffer.ini += inc_val; \
	buffer.alpha_ik += MACROVAR__nN_inc_val; \
}

#define INC_GCNMOPS_D_BUF(buffer, inc_val) \
{ \
	buffer.lambda_est+=inc_val; \
	buffer.alpha_est+=inc_val; \
	buffer.expCN_idx+=inc_val; \
	buffer.ExpLogFoldChange+=inc_val; \
	buffer.ini+=inc_val; \
	buffer.alpha_ik+=inc_val; \
}

#define INC_GCNMOPS_SXP_DATA_PTRS(inc_val, sxp_dataPtrs) \
{ \
	sxp_dataPtrs.lambda_est			+= inc_val; \
	sxp_dataPtrs.alpha_est			+= inc_val; \
	sxp_dataPtrs.expCN_idx			+= inc_val; \
	sxp_dataPtrs.ExpLogFoldChange	+= inc_val; \
	sxp_dataPtrs.ini				+= inc_val; \
	sxp_dataPtrs.alpha_ik			+= inc_val;	\
}



#endif
