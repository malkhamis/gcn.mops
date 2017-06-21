#ifndef GCNMOPS_PRIV_H
#define GCNMOPS_PRIV_H

#include "gcnmops_typedefs.h"

// returns a table of chunks in which each chunk stores pointer to a GRange submatrix
struct X_partition_table _create_partition_table(double *matrix, unsigned int x_nRanges, unsigned int nSamples, unsigned int chunk_len);
// frees internal members of table
void _free_partition_table(struct partition_table *table);


// alloc space for results buffer and set pointers
struct gcnmops_buf _alloc_gcnmops_buf(char d_or_h, unsigned int buf_len, unsigned int numClasses, unsigned int numSamples);
// frees internal members of 'buffer'
void _free_gcnmops_buf(char d_or_h, struct gcnmops_buf buffer);


// allocate space for args and constants that are common for all x_norm rows
struct gcnmops_args _alloc_gcnmops_args(char d_or_h, unsigned int numClasses, unsigned int numSamples);
// free internal members of 'args' and other precomputed args-dependent consts
void _free_gcnmops_args(char d_or_h, struct gcnmops_args args);


// pre-alloc space for gcnmops's internal tmp vars
struct gcnmops_tmpVars _prealloc_gcnmops_workspace(char d_or_h, unsigned int nTmpSpaces, unsigned int numClasses, unsigned int numSamples);
// frees internal members of 'workspace'
void _free_gcnmops_workspace(char d_or_h, struct gcnmops_tmpVars workspace);

#endif
