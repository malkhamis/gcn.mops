#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

using namespace std;


extern "C" SEXP cnmops(SEXP xS, SEXP IS, SEXP covS, SEXP cycS, SEXP alphaInitS,
		SEXP lambdaInitS, SEXP alphaPriorS) {
	double eps=1e-100;
	int i, k;
	int N = length(xS);
	int n = length(IS);
	int cyc = (int)(INTEGER(cycS)[0]);

	double *x=REAL(xS);
	double *I=REAL(IS);
	double *alphaInit=REAL(alphaInitS);
	double *lambdaInit=REAL(lambdaInitS);
	double *alphaPrior=REAL(alphaPriorS);
	double *cov=REAL(covS);

	double meanx=0.0;
	double *lg=Calloc(N, double);
	for(k = 0; k < N; k++) {
		lg[k] = lgammafn(x[k]+1);
		meanx += x[k];
	}
	meanx=meanx/N;

	double sumAlphaPrior=0.0;
	for(i = 0; i < n; i++)
		sumAlphaPrior=sumAlphaPrior+alphaPrior[i];

	SEXP alpha_ik_RET;
	PROTECT(alpha_ik_RET = allocMatrix(REALSXP, n, N));
	double *alpha_ik=REAL(alpha_ik_RET);

	SEXP alpha_i_RET;
	PROTECT(alpha_i_RET = allocVector(REALSXP, n));
	double *alpha_i=REAL(alpha_i_RET);

	SEXP alpha_est_RET;
	PROTECT(alpha_est_RET = allocVector(REALSXP, n));
	double *alpha_est=REAL(alpha_est_RET);

	SEXP lambda_est_RET;
	PROTECT(lambda_est_RET = allocVector(REALSXP, n));
	double *lambda_est=REAL(lambda_est_RET);

	for(i = 0; i < n; i++) {
		alpha_est[i]=alphaInit[i];
		lambda_est[i]=lambdaInit[i]*I[i];
	}

	double colSum=0.0;
	double sumIalpha_i=0.0;
	for(int cycNr=0; cycNr < cyc; cycNr++) {
		//k==row, N==col
                //Rprintf("alpha_ik: \n");
		for (k=0; k<N; k++) {
			colSum=0;
			for (i=0; i<n; i++){
				alpha_ik[k*n+i] = alpha_est[i]*exp((x[k]*
						log(cov[k]*lambda_est[i]))-lg[k]-cov[k]*lambda_est[i]);
				if (alpha_ik[k*n+i] < eps){
					alpha_ik[k*n+i] = eps;
				}
				colSum=colSum+alpha_ik[k*n+i];
			}
			for (i=0; i<n; i++) {
				alpha_ik[k*n+i]=alpha_ik[k*n+i]/(colSum);
				//Rprintf("%lf   ", alpha_ik[k*n+i]);
			}
			//Rprintf("\n");
		}

		sumIalpha_i=0.0;
		for (i=0; i<n; i++) {
			alpha_i[i]=0.0;
			for (k=0; k<N; k++){
				alpha_i[i]=alpha_i[i]+alpha_ik[k*n+i];
				sumIalpha_i=sumIalpha_i+I[i]*alpha_ik[k*n+i]*cov[k];
			}
			alpha_i[i]=alpha_i[i]/N;

		}

		for(i=0; i<n; i++) {
			lambda_est[i]=N*meanx/sumIalpha_i*I[i];
			alpha_est[i]=(alphaPrior[i]+alpha_i[i])/(1+sumAlphaPrior);
		}
		/*Rprintf("lambda_est: \n");
		//for(i=0; i<n; i++) {
		//  Rprintf("%lf   ", lambda_est[i]);
		//}
		//Rprintf("\n");
                Rprintf("alpha_est: \n");
                for(i=0; i<n; i++) {
                  Rprintf("%lf   ", alpha_est[i]);
		    }
                Rprintf("\n");
		*/

	}

	Free(lg);

	SEXP namesRET;
	PROTECT(namesRET = allocVector(STRSXP, 4));
	SET_STRING_ELT(namesRET, 0, mkChar("alpha.ik"));
	SET_STRING_ELT(namesRET, 1, mkChar("alpha.i"));
	SET_STRING_ELT(namesRET, 2, mkChar("alpha.est"));
	SET_STRING_ELT(namesRET, 3, mkChar("lambda.est"));

	SEXP RET;
	PROTECT(RET = allocVector(VECSXP, 4));
	SET_VECTOR_ELT(RET, 0, alpha_ik_RET);
	SET_VECTOR_ELT(RET, 1, alpha_i_RET);
	SET_VECTOR_ELT(RET, 2, alpha_est_RET);
	SET_VECTOR_ELT(RET, 3, lambda_est_RET);
	setAttrib(RET, R_NamesSymbol, namesRET);
	UNPROTECT(6);
	return(RET);
}
