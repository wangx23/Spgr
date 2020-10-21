#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _Spgr_BIC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_BICc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_BICcr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_BICcrx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_BICcx(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_BIClog(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_BIClogr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_BIClogrx(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_BIClogx(SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_BICr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_BICrx(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_BICx(SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_cal_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_cal_initial0(SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_cal_initialr(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_cal_initialr2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_cal_initialrx(SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_cal_initialrx2(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_cal_initialx(SEXP, SEXP, SEXP);
extern SEXP _Spgr_getgroup(SEXP, SEXP, SEXP);
extern SEXP _Spgr_getorder(SEXP);
extern SEXP _Spgr_ngetgroup(SEXP, SEXP, SEXP);
extern SEXP _Spgr_selectlam(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_selectlamx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_Spgr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_Spgr_penalty(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_Spgrr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_Spgrrx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_Spgrwise_lasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_Spgrwise_rep_lasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_Spgrwise_rep_scad(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_Spgrwise_scad(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Spgr_Spgrx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_Spgr_BIC",                (DL_FUNC) &_Spgr_BIC,                 5},
    {"_Spgr_BICc",               (DL_FUNC) &_Spgr_BICc,                6},
    {"_Spgr_BICcr",              (DL_FUNC) &_Spgr_BICcr,               7},
    {"_Spgr_BICcrx",             (DL_FUNC) &_Spgr_BICcrx,              6},
    {"_Spgr_BICcx",              (DL_FUNC) &_Spgr_BICcx,               5},
    {"_Spgr_BIClog",             (DL_FUNC) &_Spgr_BIClog,              5},
    {"_Spgr_BIClogr",            (DL_FUNC) &_Spgr_BIClogr,             6},
    {"_Spgr_BIClogrx",           (DL_FUNC) &_Spgr_BIClogrx,            5},
    {"_Spgr_BIClogx",            (DL_FUNC) &_Spgr_BIClogx,             4},
    {"_Spgr_BICr",               (DL_FUNC) &_Spgr_BICr,                6},
    {"_Spgr_BICrx",              (DL_FUNC) &_Spgr_BICrx,               5},
    {"_Spgr_BICx",               (DL_FUNC) &_Spgr_BICx,                4},
    {"_Spgr_cal_initial",        (DL_FUNC) &_Spgr_cal_initial,         4},
    {"_Spgr_cal_initial0",       (DL_FUNC) &_Spgr_cal_initial0,        4},
    {"_Spgr_cal_initialr",       (DL_FUNC) &_Spgr_cal_initialr,        5},
    {"_Spgr_cal_initialr2",      (DL_FUNC) &_Spgr_cal_initialr2,       6},
    {"_Spgr_cal_initialrx",      (DL_FUNC) &_Spgr_cal_initialrx,       4},
    {"_Spgr_cal_initialrx2",     (DL_FUNC) &_Spgr_cal_initialrx2,      5},
    {"_Spgr_cal_initialx",       (DL_FUNC) &_Spgr_cal_initialx,        3},
    {"_Spgr_getgroup",           (DL_FUNC) &_Spgr_getgroup,            3},
    {"_Spgr_getorder",           (DL_FUNC) &_Spgr_getorder,            1},
    {"_Spgr_ngetgroup",          (DL_FUNC) &_Spgr_ngetgroup,           3},
    {"_Spgr_selectlam",          (DL_FUNC) &_Spgr_selectlam,          11},
    {"_Spgr_selectlamx",         (DL_FUNC) &_Spgr_selectlamx,         10},
    {"_Spgr_Spgr",               (DL_FUNC) &_Spgr_Spgr,               11},
    {"_Spgr_Spgr_penalty",       (DL_FUNC) &_Spgr_Spgr_penalty,       13},
    {"_Spgr_Spgrr",              (DL_FUNC) &_Spgr_Spgrr,              12},
    {"_Spgr_Spgrrx",             (DL_FUNC) &_Spgr_Spgrrx,             11},
    {"_Spgr_Spgrwise_lasso",     (DL_FUNC) &_Spgr_Spgrwise_lasso,     10},
    {"_Spgr_Spgrwise_rep_lasso", (DL_FUNC) &_Spgr_Spgrwise_rep_lasso, 11},
    {"_Spgr_Spgrwise_rep_scad",  (DL_FUNC) &_Spgr_Spgrwise_rep_scad,  12},
    {"_Spgr_Spgrwise_scad",      (DL_FUNC) &_Spgr_Spgrwise_scad,      11},
    {"_Spgr_Spgrx",              (DL_FUNC) &_Spgr_Spgrx,              10},
    {NULL, NULL, 0}
};

void R_init_Spgr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
