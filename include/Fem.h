#ifndef __header_Fem_h__
#define __header_Fem_h__

#include "eFace.h"
#include "Workelm.h"
#include "ns_config.h"
#include "Err.h"

#define _build_partialoff 		( _nu*_nu*_nu*_ndg * 3 + _nu*_ndg*_ndg*4   +  _nu*_np   )
#define _build_twice_partialoff  	( _nu*_nu*_ndg*2+_nu*_ndg*_ndg*2)

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct
  {
    eFace elm;
 
    /* ORDER MATTERS */
    pR exact_dyu_p;
    pR exact_dens_dyu_u_u;
    pR exact_dens_u_dyu_u;
    pR exact_u_dyu_dens_u;
    pR exact_u_dyf_f;
    pR exact_dyf_u_f;
    pR exact_u_f_dyf;
    pR exact_f_u_dyf;

    pR exact_dxu_p;
    pR exact_dens_dxu_u_u;
    pR exact_dens_u_dxu_u;
    pR exact_u_dxu_dens_u;
    pR exact_u_dxf_f;
    pR exact_dxf_u_f;
    pR exact_u_f_dxf;
    pR exact_f_u_dxf;

    pR exact_visc_dxu_dxu;
    pR exact_dxu_visc_dxu;
    pR exact_u_dxf_dxf;
    pR exact_dxf_u_dxf;

    pR exact_visc_dxu_dyu;
    pR exact_dxu_visc_dyu;
    pR exact_u_dxf_dyf;
    pR exact_dxf_u_dyf;

    pR exact_visc_dyu_dxu;
    pR exact_dyu_visc_dxu;
    pR exact_u_dyf_dxf;
    pR exact_dyf_u_dxf;

    pR exact_visc_dyu_dyu;
    pR exact_dyu_visc_dyu;
    pR exact_u_dyf_dyf;
    pR exact_dyf_u_dyf;

    cst_pR exact_u_dens_u;
    cst_pR exact_dens_u_u;

    R 		J_UU		[_nu*_nu];
    R 		J_UV		[_nu*_nu];
    R 		J_VU		[_nu*_nu];
    R 		J_VV		[_nu*_nu];

    R 		visc_dxu_dxu	[_nu*_nu];
    R 		visc_dyu_dyu	[_nu*_nu];
    R 		visc_dxu_dyu	[_nu*_nu];
    R 		visc_dyu_dxu	[_nu*_nu];

    R 		dxu_visc_dxu	[_nu*_ndg];
    R 		dyu_visc_dyu	[_nu*_ndg];
    R 		dxu_visc_dyu	[_nu*_ndg];
    R 		dyu_visc_dxu	[_nu*_ndg];

    R 		dxv_visc_dxv	[_nu*_ndg];
    R 		dyv_visc_dyv	[_nu*_ndg];
    R 		dxv_visc_dyv	[_nu*_ndg];
    R 		dyv_visc_dxv	[_nu*_ndg];

    R 		dens_u_u	[_nu*_nu];
    R 		dens_dxu_u	[_nu*_nu];
    R 		dens_dyu_u	[_nu*_nu];
    R 		u_dens_u	[_ndg*_nu];
    R 		v_dens_v	[_ndg*_nu];
    R 		dens_dyu_v_u	[_ndg*_nu];
    R 		dens_v_dyu_u	[_nu*_nu*_nu];
    R 		dens_u_dyu_u	[_nu*_nu*_nu];
    R 		dens_dxu_u_u	[_nu*_nu*_nu];
    R 		dens_u_dxu_u	[_nu*_nu*_nu];
    R 		dens_dyu_u_u	[_nu*_nu*_nu];

    R 		build_partialx		[2*_build_partialoff];
    R 		build_partialr		[2*_build_partialoff];
  
    R 		build_twice_partialx	[4*_build_twice_partialoff];
    R 		build_twice_partialr	[4*_build_twice_partialoff];

    R 		shaper		[3*128];
    R 		shapeu		[_nu*128];
    R 		drshapeu	[_nu*128];
    R 		dsshapeu	[_nu*128];
    R shapep		[_np*128];
    R drshapep	[_np*128];
    R dsshapep	[_np*128];
    R shapedg	[_ndg*128];
    R drshapedg	[_ndg*128];
    R dsshapedg	[_ndg*128];

    R qref_radius_dens_u_u		[3*_nu*_nu*_ndg];

    R qref_radius_visc_dru_dru	[3*_nu*_nu*_ndg];
    R qref_radius_visc_dru_dsu	[3*_nu*_nu*_ndg];
    R qref_radius_visc_dsu_dru	[3*_nu*_nu*_ndg];
    R qref_radius_visc_dsu_dsu	[3*_nu*_nu*_ndg];

    R qref_radius_p_dsu		[3*_nu*_np];
    R qref_radius_p_dru		[3*_nu*_np];
    R qref_radius_dsu_p		[3*_nu*_np];
    R qref_radius_dru_p		[3*_nu*_np];

    R qref_dens_u_u		[3*_nu*_nu*_ndg];
    R qref_visc_dru_dru	[3*_nu*_nu*_ndg];
    R qref_visc_dru_dsu	[3*_nu*_nu*_ndg];
    R qref_visc_dsu_dru	[3*_nu*_nu*_ndg];
    R qref_visc_dsu_dsu	[3*_nu*_nu*_ndg];
    R qref_p_dsu		[3*_nu*_np];
    R qref_p_dru		[3*_nu*_np];
    R qref_dsu_p		[3*_nu*_np];
    R qref_dru_p		[3*_nu*_np];
  
  } Fem,*RESTRICT pFem;

  typedef const Fem * RESTRICT pFemReadOnly;


  cst_pR  FemReadOnly_get_dens_u_u	(pFemReadOnly	const 	self_);
  pR  	Fem_get_dens_u_u		(pFem		const 	self_);


  cst_pR  FemReadOnly_getJ_UU		(pFemReadOnly	const 	self_);
  pR  	Fem_getJ_UU			(pFem		const 	self_);

  cst_pR  FemReadOnly_getJ_VV		(pFemReadOnly	const 	self_);
  pR  	Fem_getJ_VV			(pFem		const 	self_);

  cst_pR  FemReadOnly_getJ_UV		(pFemReadOnly	const 	self_);
  pR  	Fem_getJ_UV			(pFem		const 	self_);

  cst_pR  FemReadOnly_getJ_VU		(pFemReadOnly	const 	self_);
  pR  	Fem_getJ_VU			(pFem		const 	self_);

  void 	Fem_setelm			(pFem		const 	self_,
					 pWorkelm 	const 	elm_);

  void 	Fem_def				(pFem 		const 	self_,
					 cst_eFace 		elm_,
					 STR 			errmsg_,
					 Err*			err_);

  void 	Fem_def_from_quadrature		(pFem 		const 	self_,
					 cst_pI			qn_,
					 cst_pR			qw_,
					 cst_pR			qp_,
					 cst_pI			qoff_,
					 STR 			errmsg_,
					 Err*			err_);

#ifdef __cplusplus
}
#endif

#endif

