// Copyright (C) 2016 Mohammad Alkhamis
// <malkhamis@protonmail.com>

#include <stdlib.h>
#include <Rdefines.h>
#include "macros.h"
#include "gcnmops_typedefs.h"
#include <cuda_runtime.h>

using namespace std;

struct X_partition_table _create_partition_table(double *matrix, unsigned int x_nRanges, unsigned int nSamples, unsigned int chunk_len) {

#ifdef DEBUG
	if(matrix == NULL)
		error("%s:%s @ line %d\nfailed to create partition table for matrix.\naddress of matrix (=%p) is invalid\n",
			__FILE__, __func__, __LINE__, matrix);
		
	if(x_nRanges < 1 || nSamples < 1 || chunk_len < 1)
		error("%s:%s @ line %d\nfailed to create partition table for X.\nnumber of GRanges (=%d), samples (=%d), and length of chunk (=%d) must be integers greater than zero\n", 
			__FILE__, __func__, __LINE__, x_nRanges, nSamples, chunk_len);
	
	if( chunk_len > x_nRanges )
		error("%s:%s @ line %d\nfailed to create partition table for matrix.\nchunk length (=%d) cannot be greater than number of GRanges (=%d) in matrix\n",
			 __FILE__, __func__, __LINE__, chunk_len, x_nRanges);
#else
	if(matrix == NULL || x_nRanges < 1 || nSamples < 1 || chunk_len < 1 || chunk_len > x_nRanges)
		error("failed to create partition table. Invalid arguments!'\n");
#endif

	struct X_partition_table table;
		
	unsigned int last_chunk_len = x_nRanges % chunk_len;
	unsigned int nChunks = (x_nRanges / chunk_len) + (last_chunk_len == 0?0:1);

	table.chunk = Calloc(nChunks, struct X_norm);
	table.nPartitions = nChunks;
	
	if(nChunks == 1) {
		table.chunk->nSamples = nSamples;
		table.chunk->nRanges = chunk_len;
		table.chunk->data = matrix;
		
	} else {	
		for(unsigned int i = 0; i < nChunks; i++) {
			(table.chunk[i]).nSamples = nSamples;
			(table.chunk[i]).nRanges = chunk_len;
			(table.chunk[i]).data = matrix + (i*chunk_len);
		}
		// make correction if last chunk has a different length
		if(last_chunk_len > 0)
			(table.chunk[nChunks-1]).nRanges = last_chunk_len;
	}
	
	return table;
}

void _free_partition_table(struct X_partition_table *table){

	Free(table->chunk);
}


struct gcnmops_buf _alloc_gcnmops_buf(char d_or_h, unsigned int buf_len, unsigned int numClasses, unsigned int numSamples) {

#ifdef DEBUG
	if(buf_len < 1 || numClasses < 1 || numSamples < 1)
		error("%s:%s @ line %d\nfailed to allocate buffer\nlength of buffer (=%d), N (=%d), and n (=%d) must be integers greater than zero\n",
			__FILE__, __func__, __LINE__, buf_len, numSamples, numClasses);

	if(!((d_or_h == 'd') || (d_or_h == 'D') || (d_or_h == 'h') || (d_or_h == 'H')))
		error("%s:%s @ line %d\nfailed to allocate buffer\ntype of memory (=%c) must either be device ('d'/'D') or host ('h'/'H')\n",
			__FILE__, __func__, __LINE__, d_or_h);	
#else
	if(buf_len < 1 || numClasses < 1 || numSamples < 1 || 
	!((d_or_h == 'd') || (d_or_h == 'D') || (d_or_h == 'h') || (d_or_h == 'H')) )
		error("failed to allocate buffer. Invalid arguments!\n");
#endif
	
	struct gcnmops_buf buffer;
	unsigned int slice_nElements = N_ELMS_GCNMOPS_BUFSLICE(numSamples, numClasses);
	unsigned long long nElements = buf_len * slice_nElements;

	double *data;
	size_t nBytes = nElements * sizeof(double);
	
	switch(d_or_h) {
		case 'h':
		case 'H':
		{
			//data = Calloc(nElements, double);
			// we instead allocate page-locked/pinned memory (faster)
			CHK_CU(cudaMallocHost((double**)&data, nBytes));
			break;
		}
		
		case 'd':
		case 'D':
		{
			CHK_CU(cudaMalloc((double**)&data, nBytes));
			break;
		}
		
		default:
			{break;}
	}
	
	// metadata
	buffer.nSlices = buf_len;
	buffer.nBytes_slice = slice_nElements * sizeof(double);
	buffer.mem_handle = data; // to free mem
	
	// wire memory space to the buffer internal pointers
	unsigned long long offset_alpha_est = 0; // length per x_norm row  = numClasses
	unsigned long long offset_lambda_est = offset_alpha_est + (numClasses*buf_len); // len = numClasses
	unsigned long long offset_alpha_ik = offset_lambda_est + (numClasses*buf_len); // length per x_norm row = numClasses*numSamples
	unsigned long long offset_ini = offset_alpha_ik + (numClasses*numSamples*buf_len); // length per x_norm row = 1
	unsigned long long offset_ExpLog = offset_ini + (1*buf_len); // length per x_norm row = numSamples
	unsigned long long offset_expCN_idx = offset_ExpLog + (numSamples*buf_len); // length per x_norm row = numSamples
	
	buffer.alpha_est		= data + offset_alpha_est;
	buffer.lambda_est		= data + offset_lambda_est;
	buffer.alpha_ik			= data + offset_alpha_ik;
	buffer.ini				= data + offset_ini;
	buffer.ExpLogFoldChange	= data + offset_ExpLog;
	buffer.expCN_idx		= data + offset_expCN_idx;
		
	return buffer;
}

void _free_gcnmops_buf(char d_or_h, struct gcnmops_buf buffer) {
	
#ifdef DEBUG
	if(!((d_or_h == 'd') || (d_or_h == 'D') || (d_or_h == 'h') || (d_or_h == 'H')))
		error("%s:%s @ line %d\nfailed to free buffer\ntype of memory (=%c) must either be device ('d'/'D') or host ('h'/'H')\n",
			__FILE__, __func__, __LINE__, d_or_h);
#else
	if(!((d_or_h == 'd') || (d_or_h == 'D') || (d_or_h == 'h') || (d_or_h == 'H')))
		error("failed to free buffer. Invalid arguments!\n");
#endif

	switch(d_or_h) {
		case 'h':
		case 'H':
		{
			//Free(buffer.mem_handle);
			CHK_CU( cudaFreeHost(buffer.mem_handle) );
			break;
		}
		
		case 'd':
		case 'D':
		{
			CHK_CU( cudaFree(buffer.mem_handle) );
			break;
		}
			
		default:
			{break;}
	}
}


struct gcnmops_args _alloc_gcnmops_args(char d_or_h, unsigned int numClasses, unsigned int numSamples) {

#ifdef DEBUG
	if(numClasses < 1)
		error("%s:%s @ line %d\nfailed to allocate memory for gcnmop's pre-copmuted constants\nvalue of n (=%d) must be unsigned integers greater than zero\n",__FILE__, __func__, __LINE__, numClasses);

	if(!((d_or_h == 'd') || (d_or_h == 'D') || (d_or_h == 'h') || (d_or_h == 'H')))
		error("%s:%s @ line %d\nfailed to allocate memory for gcnmop's pre-copmuted constants\ntype of memory (=%c) must either be device ('d'/'D') or host ('h'/'H')\n", __FILE__, __func__, __LINE__, d_or_h);	
#else
	if(numClasses < 1 || 
	!((d_or_h == 'd') || (d_or_h == 'D') || (d_or_h == 'h') || (d_or_h == 'H')) )
		error("failed to allocate memory for gcnmop's pre-copmuted constants. Invalid arguments!\n");
#endif
	
	unsigned int nElements = N_ELMS_GCNMOPS_ARR_ARGS(numSamples, numClasses);
	double *linear_space;
	
	switch(d_or_h) {
		case 'h': // in case we need for one reason or another
		case 'H':
		{
			linear_space = Calloc(nElements, double);
			break;
		}
		
		case 'd':
		case 'D':
		{
			size_t nBytes = nElements * sizeof(double);
			CHK_CU( cudaMalloc((double**)&linear_space, nBytes) );
			break;
		}
		default:
			{break;}
	}
	
	// wire the internal struct members to the mem space
	struct gcnmops_args args;
	
	args.mem_handle = linear_space; // to free mem
	
	args.I = linear_space;
	args.cov = args.I + numClasses;
	args.alphaInit = args.cov + numSamples;
	args.alphaPrior = args.alphaInit + numClasses;
	
	args.log2_I = args.alphaPrior + numClasses;
	args.abs_log2_I = args.log2_I + numClasses;

	return args;
}


void _free_gcnmops_args(char d_or_h, struct gcnmops_args args) {

#ifdef DEBUG
	if(!((d_or_h == 'd') || (d_or_h == 'D') || (d_or_h == 'h') || (d_or_h == 'H')))
		error("%s:%s @ line %d\nfailed to free constants' memory\ntype of memory (=%c) must either be device ('d'/'D') or host ('h'/'H')\n",
			__FILE__, __func__, __LINE__, d_or_h);
#else
	if(!((d_or_h == 'd') || (d_or_h == 'D') || (d_or_h == 'h') || (d_or_h == 'H')))
		error("failed to free buffer. Invalid arguments!\n");
#endif

	switch(d_or_h) {
		case 'h':
		case 'H':
		{
			Free(args.mem_handle);
			break;
		}
		
		case 'd':
		case 'D':
		{
			CHK_CU( cudaFree(args.mem_handle) );
			break;
		}
		
		default:
			{break;}
	}
	
}


// each gpu thread will have a pre-allocated mem space for tmp vars
struct gcnmops_tmpVars _prealloc_gcnmops_workspace(char d_or_h, unsigned int nTmpSpaces, unsigned int numClasses, unsigned int numSamples) {

#ifdef DEBUG
	if(nTmpSpaces < 1 || numClasses < 1 || numSamples < 1)
		error("%s:%s @ line %d\nfailed to pre-allocate memory for gcnmop's internal temporary variables\nnumber of workspaces (=%d), N (=%d), and n (=%d) must be unsigned integers greater than zero\n",
			__FILE__, __func__, __LINE__, nTmpSpaces, numSamples, numClasses);

	if(!((d_or_h == 'd') || (d_or_h == 'D') || (d_or_h == 'h') || (d_or_h == 'H')))
		error("%s:%s @ line %d\nfailed to pre-allocate memory for gcnmop's internal temporary variables\ntype of memory (=%c) must either be device ('d'/'D') or host ('h'/'H')\n",
			__FILE__, __func__, __LINE__, d_or_h);	
#else
	if(nTmpSpaces < 1 || numClasses < 1 || numSamples < 1 || 
	!((d_or_h == 'd') || (d_or_h == 'D') || (d_or_h == 'h') || (d_or_h == 'H')) )
		error("failed to pre-allocate memory for gcnmop's internal temporary variables. Invalid arguments!\n");
#endif

	unsigned long long nElements = nTmpSpaces * N_ELMS_GCNMOPS_TMPVARS(numSamples, numClasses);
	double *linear_space;
	
	switch(d_or_h) {
		case 'h':
		case 'H':
		{
			linear_space = Calloc(nElements, double);
			break;
		}
		
		case 'd':
		case 'D':
		{
			unsigned long long nBytes_elements = nElements * sizeof(double);
			CHK_CU( cudaMalloc((double**)&linear_space, nBytes_elements) );
			break;
		}
		
		default:
			{break;}
	}

	// offsets to wire every workspace to its memory space
	// current layout: ===x_over_cov_sorted===lambdaInit===lg===alpha_i===draftspace_ini
	unsigned long long offset_x_over_cov_sorted = 0; //in case we change the layout
	unsigned long long offset_lambdaInit = offset_x_over_cov_sorted + (numSamples * nTmpSpaces);
	unsigned long long offset_lg = offset_lambdaInit + (numClasses * nTmpSpaces);
	unsigned long long offset_alpha_i = offset_lg + (numSamples * nTmpSpaces);
	unsigned long long offset_draftspace_ini = offset_alpha_i + (numClasses * nTmpSpaces);
		
	struct gcnmops_tmpVars workspace;
	workspace.mem_handle = linear_space;
	
	workspace.x_over_cov_sorted = linear_space + offset_x_over_cov_sorted;
	workspace.lambdaInit = linear_space + offset_lambdaInit;
	workspace.lg = linear_space + offset_lg;
	workspace.alpha_i = linear_space + offset_alpha_i;
	workspace.draftspace_ini = linear_space + offset_draftspace_ini;
	
	return(workspace);
}


void _free_gcnmops_workspace(char d_or_h, struct gcnmops_tmpVars workspace) {
	
#ifdef DEBUG
	if(!((d_or_h == 'd') || (d_or_h == 'D') || (d_or_h == 'h') || (d_or_h == 'H')))
		error("%s:%s @ line %d\nfailed to free buffer\ntype of memory (=%c) must either be device ('d'/'D') or host ('h'/'H')\n",
			__FILE__, __func__, __LINE__, d_or_h);
#else
	if(!((d_or_h == 'd') || (d_or_h == 'D') || (d_or_h == 'h') || (d_or_h == 'H')))
		error("failed to free buffer. Invalid arguments!\n");
#endif

	switch(d_or_h) {
		case 'h':
		case 'H':
		{
			Free(workspace.mem_handle);
			break;
		}
		
		case 'd':
		case 'D':
		{
			CHK_CU( cudaFree(workspace.mem_handle) );
			break;
		}
		
		default:
			{break;}
	}
}
