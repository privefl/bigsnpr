#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP bigsnpr_clumping(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bigsnpr_clumping2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bigsnpr_corMat(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bigsnpr_doubleBM(SEXP, SEXP);
extern SEXP bigsnpr_local_clumping(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bigsnpr_pruning(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bigsnpr_pruning2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bigsnpr_readbina(SEXP, SEXP, SEXP);
extern SEXP bigsnpr_roll_mean(SEXP, SEXP);
extern SEXP bigsnpr_testWrite(SEXP, SEXP);
extern SEXP bigsnpr_writebina(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"bigsnpr_clumping",       (DL_FUNC) &bigsnpr_clumping,        9},
  {"bigsnpr_clumping2",      (DL_FUNC) &bigsnpr_clumping2,      10},
  {"bigsnpr_corMat",         (DL_FUNC) &bigsnpr_corMat,          5},
  {"bigsnpr_doubleBM",       (DL_FUNC) &bigsnpr_doubleBM,        2},
  {"bigsnpr_local_clumping", (DL_FUNC) &bigsnpr_local_clumping,  7},
  {"bigsnpr_pruning",        (DL_FUNC) &bigsnpr_pruning,         9},
  {"bigsnpr_pruning2",       (DL_FUNC) &bigsnpr_pruning2,       10},
  {"bigsnpr_readbina",       (DL_FUNC) &bigsnpr_readbina,        3},
  {"bigsnpr_roll_mean",      (DL_FUNC) &bigsnpr_roll_mean,       2},
  {"bigsnpr_testWrite",      (DL_FUNC) &bigsnpr_testWrite,       2},
  {"bigsnpr_writebina",      (DL_FUNC) &bigsnpr_writebina,       3},
  {NULL, NULL, 0}
};

void R_init_bigsnpr(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
