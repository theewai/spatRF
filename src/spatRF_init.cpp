#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


extern "C"{
extern SEXP C_dist(SEXP, SEXP, SEXP);
extern SEXP C_splits(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_splits1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"C_dist",  (DL_FUNC) &C_dist,  3},
  {"C_splits",  (DL_FUNC) &C_splits,  12},
  {"C_splits1", (DL_FUNC) &C_splits1,  9},
  {NULL, NULL, 0}
};

void R_init_spatRF(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
}