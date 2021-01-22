#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
  /* .C calls */
  extern void gutsredit_free();
extern void gutsredsd_free();

static const R_CMethodDef CEntries[] = {
  {"gutsredit_free", (DL_FUNC) &gutsredit_free, 0},
  {"gutsredsd_free", (DL_FUNC) &gutsredsd_free, 0},
  {NULL, NULL, 0}
};

void R_init_morse(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}