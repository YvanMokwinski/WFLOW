#ifndef __header_Workelm_h__
#define __header_Workelm_h__
#include "eDim.h"
#include "ns_config.h"
#include "ensBASIS.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct
  {
    R 		locmat	 [128*128];
    R 		locresidu[128*128];
    I 		locnumer [128];
    R 		alpha2;
    R 		ialpha2;
    R  		edgelm[__eDim_ALL*8];
    R 		normalelmx[__eDim_ALL*8];
    R 		normalelmy[__eDim_ALL*8];
    R 		alpha3[__eDim_ALL+1];
    R		belm	[__eDim_ALL*__eDim_ALL];
    R		sbelm	[__eDim_ALL*__eDim_ALL];
    R 		cooelm	[16];
    R 		dgelm[_ndg];
    R 		dgelmi[_ndg];
    R 		dgelmii[_ndg];
    cst_pR		uelm;
    cst_pR		velm;
    cst_pR		welm;
    cst_pR		pelm;
    cst_pR		uelmi;
    cst_pR		velmi;
    cst_pR		welmi;
    cst_pR		pelmi;
    cst_pR		uelmii;
    cst_pR		velmii;
    cst_pR		welmii;
    cst_pR		pelmii;
#define __NS_VAR_NMAX__ 4
    R 		denselm		[__NS_SHAPE_NMAX__];
    R 		viscelm		[__NS_SHAPE_NMAX__];
    R 		denselmi	[__NS_SHAPE_NMAX__];
    R 		viscelmi	[__NS_SHAPE_NMAX__];
    R 		kappaelm	[__NS_SHAPE_NMAX__];
    R 		ddlelm		[__NS_SHAPE_NMAX__*__NS_VAR_NMAX__];
    R 		ddlelmi		[__NS_SHAPE_NMAX__*__NS_VAR_NMAX__];
    R 		ddlelmii	[__NS_SHAPE_NMAX__*__NS_VAR_NMAX__];
  }  Workelm,*RESTRICT pWorkelm;

  typedef const Workelm * pWorkelmReadOnly;

  void Workelm_init(pWorkelm const self_);
  
#ifdef __cplusplus
}
#endif

#endif
