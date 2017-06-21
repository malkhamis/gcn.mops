#include <Rdefines.h>

/* cnmops_GPU_wrapper.cpp */
SEXP gcnmops_w(SEXP gpuS, SEXP x_normS, SEXP IS, SEXP classesS, SEXP covS, SEXP cycS, SEXP idxCN2S, SEXP alphaInitS,
		SEXP alphaPriorS, SEXP minReadCountS, SEXP returnPosteriorS);

/* segment.cpp */
SEXP segment(SEXP xS, SEXP epsS, SEXP deltaS, SEXP maxIntS, SEXP minSegS, SEXP squashingS, SEXP cyberWeightS);

