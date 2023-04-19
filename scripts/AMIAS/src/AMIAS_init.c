// Automatically generated, editing not advised.
#ifndef R_AMIAS_H
#define R_AMIAS_H
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("AMIAS", String)
#else
#define _(String) (String)
#endif

#define FDEF(name)  {#name, (DL_FUNC) &F77_SUB(name), sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
void F77_SUB(fusedl0)(
            int *n,
			double *y,	
			double *beta,
			double *z,
			double *u,
      int *A,
      int *I,
      int *Anull,
			int *T,
			double *rho,
			int *iter,
			int *itermax,
			int *adjust,
			int *adjust_max,
			int *delta 
			);
static R_NativePrimitiveArgType fusedl0_t[] = {
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(afusedl0)(
            int *n,
			double *y,	
			double *beta,
			double *z,
			double *u,
			int *tao,
			int *Kmax,
			int *L,
			double *eps,
			double *rho,
			int *miter,
			int *itermax,
			int *adjust,
			int *adjust_max,
			int *delta 
			);
static R_NativePrimitiveArgType afusedl0_t[] = {
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(btfusedl0)(
            int *n,
			double *y,	
			double *beta,
			double *z,
			double *u,
      int *A,
      int *I,
      int *Anull,
			int *T,
			double *rho,
			int *iter,
			int *itermax,
			double *inv,
			double *vec,
			int *veck,
			int *adjust,
			int *adjust_max,
			int *delta 
			);
static R_NativePrimitiveArgType btfusedl0_t[] = {
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(batfusedl0)(
            int *n,
			double *y,	
			double *beta,
			double *z,
			double *u,
			int *tao,
			int *Kmax,
			int *L,
			double *eps,
			double *rho,
			int *miter,
			int *itermax,
			double *inv,
			double *vec,
			int *veck,
			int *adjust,
			int *adjust_max,
			int *delta 
			);
static R_NativePrimitiveArgType batfusedl0_t[] = {
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(gfusedl0)(
            int *n,
			double *y,	
			double *beta,
			double *z,
			double *u,
      int *A,
      int *I,
      int *Anull,
			int *T,
			double *rho,
			int *iter,
			int *itermax,
			double *inv,
			double *d,
			int *m,
			int *adjust,
			int *adjust_max,
			int *delta 
			);
static R_NativePrimitiveArgType gfusedl0_t[] = {
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(agfusedl0)(
            int *n,
			double *y,	
			double *beta,
			double *z,
			double *u,
			int *tao,
			int *Kmax,
			int *L,
			double *eps,
			double *rho,
			int *miter,
			int *itermax,
			double *inv,
			double *d,
			int *m,
			int *adjust,
			int *adjust_max,
			int *delta 
			);
static R_NativePrimitiveArgType agfusedl0_t[] = {
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(tfusedl02d)(
            int *dim1,
			int *dim2,
			double *y,
			double *beta,
			double *z,
			double *u,
			int *T,
			double *rho,
			int *iter,
			int *itermax,
			double *inv,
			double *vec,
			int *veck
			);
static R_NativePrimitiveArgType tfusedl02d_t[] = {
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP
};
void F77_SUB(atfusedl02d)(
            int *dim1,
			int *dim2,
			double *y,	
			double *beta,
			double *z,
			double *u,
			int *tao,
			int *Kmax,
			int *L,
			double *eps,
			double *rho,
			int *miter,
			int *itermax,
			double *inv,
			double *vec,
			int *veck
			);
static R_NativePrimitiveArgType atfusedl02d_t[] = {
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP
};
void F77_SUB(comggfusedl0)(
            int *n,
			double *y,
			double *beta,
			double *alpha,
			double *u,
			int *T,
			double *rho,
			int *iter,
			int *itermax,
			double *inv,
			double *gamma,
			double *v,
			double *d,
			int *m,
			double *invw,
			double *W,
			int *mw,
			int *Tw,
			double *rhow
			);
static R_NativePrimitiveArgType comggfusedl0_t[] = {
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP
};
void F77_SUB(acomggfusedl0)(
            int *n,
			double *y,
			double *beta,
			double *alpha,
			double *u,
			int *tao,
			int *Kmax,
			int *L,
			double *eps,
			double *rho,
			int *miter,
			int *itermax,
			double *gamma,
			double *v,
			double *inv,
			double *d,
			int *m,
			double *invw,
			double *W,
			int *mw,
			int *Tw,
			double *rhow
			);
static R_NativePrimitiveArgType acomggfusedl0_t[] = {
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP
};
void F77_SUB(comdifusedl0)(
            int *n,
			double *y,
			double *beta,
			double *alpha,
			double *u,
			int *T,
			double *rho,
			int *iter,
			int *itermax,
			double *inv,
			double *gamma,
			double *v,
			int *Tw,
			double *rhow
			);
static R_NativePrimitiveArgType comdifusedl0_t[] = {
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP
};
void F77_SUB(acomdifusedl0)(
            int *n,
			double *y,	
			double *beta,
			double *alpha,
			double *u,
			int *tao,
			int *Kmax,
			int *L,
			double *eps,
			double *rho,
			int *miter,
			int *itermax,
			double *gamma,
			double *v,
			int *Tw,
			double *rhow
			);
static R_NativePrimitiveArgType acomdifusedl0_t[] = {
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP
};
void F77_SUB(dtdmul)(
			double *vec,
            int *n,
			int *k,
			double *mat
);
static R_NativePrimitiveArgType dtdmul_t[] = {
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP
};

static R_FortranMethodDef fMethods[] = {
  FDEF(fusedl0),
  FDEF(afusedl0),
  FDEF(btfusedl0),
  FDEF(batfusedl0),
  FDEF(gfusedl0),
  FDEF(agfusedl0),
  FDEF(tfusedl02d),
  FDEF(atfusedl02d),
  FDEF(comggfusedl0),
  FDEF(acomggfusedl0),
  FDEF(comdifusedl0),
  FDEF(acomdifusedl0),
  FDEF(dtdmul),
  {NULL, NULL, 0}
};

void R_init_AMIAS(DllInfo *dll){
  R_registerRoutines(dll, NULL, NULL, fMethods, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
#endif
