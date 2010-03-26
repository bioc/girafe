// initialization of the package

#include "girafe.h"

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/RConverters.h>
#include <R_ext/Rdynload.h>
//#include <R_ext/Utils.h>

/* LOGISTICS MOSTLY IMPORTANT FOR WINDOWS DLL */

static R_CallMethodDef girafe_calls[] = {
  {"_coverage", (DL_FUNC) &girafe_coverage, 4},
  // necessary last entry of R_CallMethodDef:
  {NULL, NULL, 0}
};

// register functions 
void R_init_girafe(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, girafe_calls, NULL, NULL);
}

void R_unload_girafe(DllInfo *dll) 
{
  // at the moment nothing to do here 
}
