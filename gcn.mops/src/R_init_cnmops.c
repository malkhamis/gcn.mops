#include "cn.mops.h"
#include "gcn.mops.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* cnmops.cpp */
	CALLMETHOD_DEF(cnmops, 7),

/* gcnmops_w.cpp */
	CALLMETHOD_DEF(gcnmops_w, 11),
	
/* segment.c */
	CALLMETHOD_DEF(segment, 7),

	{NULL, NULL, 0}
};


void R_init_cnmops(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	return;
}
