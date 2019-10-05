#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* declarations to register native routines in this package */ 

/* .C calls */
extern void wdCOV(void *, void *, void *, void *, void *, void *,void *, void *);
extern void wdCOVtest(void *, void *, void *, void *, void *, void *, void *,void *, void *,void *);

static const R_CMethodDef CEntries[] = {
  {"wdCOV",        (DL_FUNC) &wdCOV,        8},
  {"wdCOVtest",    (DL_FUNC) &wdCOVtest,    10},
  {NULL, NULL, 0}
};

void R_init_WdCor(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries,NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
