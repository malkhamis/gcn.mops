#ifndef GCNMOPS_TYPEDEFS_H
#define GCNMOPS_TYPEDEFS_H

extern "C" struct X_norm {
	unsigned int nSamples;
	unsigned int nRanges;
	double *data;
};
extern "C" struct X_partition_table {
	unsigned int nPartitions;
	struct X_norm *chunk;
};

// results will be stored contiguously in a buffer
extern "C" struct gcnmops_buf {
	// metadata
	unsigned int nSlices;
	unsigned int nBytes_slice;
	double *mem_handle; // internal use
	
	double	*alpha_est;
	double	*lambda_est;
	double	*alpha_ik;
	double	*ini;
	double	*ExpLogFoldChange;
	double	*expCN_idx;
};

// constants applicable to processing all x_nrom rows
extern "C" struct gcnmops_args {
	double *mem_handle; // internal use

	unsigned int nRanges_tot;
	unsigned int nSamples;
	unsigned int nClasses;
	//unsigned int minReadCount;
	unsigned int cyc;
	unsigned int idxCN2;
	double minReadCount;
	// arrays
	double *I;
	double *cov;
	double *alphaInit;
	double *alphaPrior;
	double *log2_I;
	double *abs_log2_I;
};

// gcnmops core's internal dynamically-allocated vars. Each GPU thread will have his tmpVars space
extern "C" struct gcnmops_tmpVars {
	double *mem_handle; //internal use
	double *x_over_cov_sorted;
	double *lambdaInit;
	double *lg;
	double *alpha_i;
	double *draftspace_ini;
};

// pointers to data pointers of SEXP objects
extern "C" struct gcnmops_SXP_DATA {
	double	*lambda_est;
	double	*alpha_est;
//	double	*expCN_idx; //unsused
	double	*ExpLogFoldChange;
	double	*ini;
	double	*alpha_ik;
};

#endif
