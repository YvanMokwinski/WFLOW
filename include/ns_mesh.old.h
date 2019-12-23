#ifndef __header_ns_mesh_h__
#define __header_ns_mesh_h__

#if 0
#include "ns_defines.h"
#endif
#include "Err.h"
#include "eDim.h"
#include "eElement.h"
#include "ensBASIS.h"


#ifdef __cplusplus
extern "C"
{
#endif


  typedef struct
  {
    eElement elm;
    /*\brief connectivity*/
    cst_pI   cnc;
    /*\brief vertex cod */
    cst_pI   cod;
    /*\brief begin vertex patchelm */
    cst_pI   bpatch;
    /*\brief index vertex patchelm */
    cst_pI   patch;
    /*\brief element adjacencies */
    cst_pI  adj;
    /*\brief coordinates */
    cst_pR  coo;
    /*\brief volume element */
    cst_pR jacelm;
    /*\brief volume edges */
    cst_pR  jacedge;
    /*\brief normaledge 			*/
    cst_pR  normaledge;
    /*\brief jacobian transformation elm */
    cst_pR trelm;
    /*\brief temps */

    /*\brief connectivity*/
    pI   own_cnc;
    /*\brief vertex cod */
    pI   own_cod;
    /*\brief begin vertex patchelm */
    pI   own_bpatch;
    /*\brief index vertex patchelm */
    pI   own_patch;
    /*\brief element adjacencies */
    pI  own_adj;
    /*\brief coordinates */
    pR  own_coo;
    /*\brief volume element */
    pR own_jacelm;
    /*\brief volume edges */
    pR  own_jacedge;
    /*\brief normaledge 			*/
    pR  own_normaledge;
    /*\brief jacobian transformation elm */
    pR own_trelm;

    pR 	time_x;  
    /*\brief nvertex */
    I nvertex;
    /*\brief nelm */
    I nelm;
    /*\brief nedge */
    I nedge;  
    /*\brief nedge boundary */
    I nedge_boundary;  
    /*\brief ntime */
    I time_n;  
    /*\brief nddl quadratic */
    /*\brief nddl quadratic bulle */
    I 	nddlP6;
    I 	noboundary_vcod;  
    eDim 	dimension;  

    void * 	spaces[__ensBASIS_n];
    R 	time_t,time_ti,time_tii;
    void * 	var_scalar	[64];
    void * 	var_vector	[64];
    void * 	var_tensor	[64];
    void * 	metrics		[64];
    I 	nddl_lagr[__ensBASIS_n];
#if 0
    nsCUBA_TRIA_ST * 		instances_cubature_triangle[16];
    nsCUBA_TRIA_INTERVAL_ST * 	instances_cubature_triangle_interval[16];
#endif
  } ns_mesh;


  
  typedef ns_mesh nsMESHZONE_ST;
  typedef nsMESHZONE_ST * nsMESHZONE;
  typedef const nsMESHZONE_ST cst_nsMESHZONE_ST;
  typedef cst_nsMESHZONE_ST* cst_nsMESHZONE;

#if 0
  void nsSPARSE_STOKES_dirichlet_init	(nsSPARSE_STOKES_ST*const 	S_,
					 const ns_mesh *const	mesh_,
					 const I 		nshape_face_u_,
					 const I 		nshape_face_p_,
					 const I 		dec_ddlu_,
					 const I 		dec_ddlv_,
					 const I 		dec_ddlp_,
					 const I 		dec_ddls_);
#endif
  
#ifndef NDEBUG
  eElement	ns_mesh_elm		(const ns_mesh*		mesh_);
  I 	ns_mesh_nelm		(const ns_mesh*		mesh_);
  I 	ns_mesh_nvertex		(const ns_mesh*		mesh_);
  I 	ns_mesh_nedge		(const ns_mesh*		mesh_);
  I 	ns_mesh_nedge_boundary	(const ns_mesh*		mesh_);
  cst_pR 	ns_mesh_trelm		(const ns_mesh*		mesh_);
  cst_pR 	ns_mesh_jacelm		(const ns_mesh*		mesh_);
  L 	ns_mesh_istransient	(const ns_mesh*		mesh_);
  eDim	ns_mesh_dimension	(const ns_mesh*		mesh_);
  I 	ns_mesh_nddlspace	(const ns_mesh* 	mesh_,
						 cst_ensBASIS 	shape_);
#else
#define 	ns_mesh_elm(_mesh) 		(_mesh)->elm
#define 	ns_mesh_nddlspace(_mesh,_shape) (_mesh)->nddl_lagr[_shape]
#define 	ns_mesh_dimension(_mesh)	(_mesh)->dimension
#define 	ns_mesh_istransient(_mesh)	( ((_mesh)->time_x) ? __ens_yes : __ens_no )
#define 	ns_mesh_nelm(_mesh) 		(_mesh)->nelm
#define 	ns_mesh_nvertex(_mesh) 		(_mesh)->nvertex
#define 	ns_mesh_nedge(_mesh) 		(_mesh)->nedge
#define 	ns_mesh_nedge_boundary(_mesh) 	(_mesh)->nedge_boundary
#define 	ns_mesh_jacelm(_mesh) 		(_mesh)->jacelm
#define 	ns_mesh_trelm(_mesh) 		(_mesh)->trelm
#endif

  void 		ns_mesh_get_trelm	(const ns_mesh*mesh_,cst_pI jelm_,pR trelm_);
  void 		ns_mesh_get_jacelm	(const ns_mesh*mesh_,cst_pI jelm_,pR jacelm_);

  void 		ns_mesh_edgeinfo		(const ns_mesh * 	mesh_,
						 cst_pI 		iedge_,
						 pI 			nadj_,
						 pI 			ielm_,
						 pI           		ilocedge_,
						 pI 			jelm_,
						 pI           		jlocedge_);

  void  ns_mesh_get_coovertex(const ns_mesh*mesh_,cst_pI ivertex_,pR xc);


  void 		ns_mesh_free			(ns_mesh * mesh_);

  void* 		ns_mesh_localization		(const ns_mesh*		mesh_);
  void* 		ns_mesh_localization_kill	(void * 		localization_);


  void 		ns_mesh_localize_recover	(void   * 		localize,
						 cst_pR		cooelm_,
						 pI 			select_,
						 pI			select_n_);

  void  		ns_mesh_localize		(void   * 		localize,
							 const ns_mesh * 	mesh_,
							 pR 			vert,
							 pI 			select_,
							 pI			select_n_);


  Err 	ns_mesh_read			(ns_mesh * 		mesh_,
						 const char * 		name_,
						 ...);


  void 		ns_mesh_cooelm			(const ns_mesh*		mesh_,
						 cst_pI 		ielm_,
						 pR         		cooelm_);

  void 		ns_mesh_patch_cooelm		(const ns_mesh*	mesh_,
						 cst_pI 	nelm_,
						 cst_pI        elm_,
						 pR         	cooelms_);


  void 		ns_mesh_cncelm			(const ns_mesh*	mesh_,
						 cst_pI 	ielm_,
						 pI         	cncelm_);

  void 		ns_mesh_cncelm0			(const ns_mesh*	mesh_,
						 cst_pI 	ielm_,
						 pI         	cncelm_);


  void 		ns_mesh_edge			(const ns_mesh * mesh_,
						 cst_pI 	iedge_,
						 pI 		nadj_,
						 pI 		ielm_,
						 pI           ilocedge_,
						 pI 		jelm_,
						 pI           jlocedge_);

  void 		ns_mesh_trelms			(const ns_mesh*	mesh_,
						 cst_pI 	start_ielm_,
						 cst_pI 	stop_ielm_,
						 pR         	trelm_,
						 cst_pI        trelmoff_);

  void 		ns_mesh_set_periodic		(ns_mesh*mesh_);


  void 		ns_mesh_extrude_scalar		(const ns_mesh * mesh_,
						 const char * 	name_,
						 const R 	coeff_,
						 const I	nframe__,
						 const I	nz_,
						 cst_pR 	z_,
						 cst_pR 	*scalars_);

  void 		ns_mesh_extrude_vector		(const ns_mesh * 	mesh_,
						 const char * 		name_,
						 const R 		coeff_,
						 const I 		nframe__,
						 const I 		nz_,
						 cst_pR 		z_,
						 cst_pR 		*scalarsx_,
						 cst_pR 		*scalarsy_);

  void 		ns_mesh_vertex_patchelm		(const ns_mesh * 	mesh_,
						 cst_pI 		ivertex_,
						 pI 			patchelm_n_,
						 pI 			patchelm_);

  void 		print_mesh_extrude_axiP2	(const ns_mesh * 	mesh_,
						 const char * 		basename_,
						 cst_pI 		nrot_,
						 cst_pI 		axe_,
						 pI 			perm_,
						 pI                    rperm_,
						 cst_pI 		nvertex_axi_,
						 cst_pI                nvertex_noaxi_);


  Err 	ns_mesh_print			(const ns_mesh * 	mesh_,
						 const char * 	name_);
  void 		ns_mesh_write_medit		(const ns_mesh*s_,
						 const char * name_,...);
  void 		ns_mesh_extrude_medit		(const ns_mesh*	mesh_,
						 const char * 	name_,
						 const I 	nframe_,
						 const I 	nz_,
						 cst_pR 	z_,
						 const R 	coeff_);

  void 		print_mesh_extrude_axi		(const ns_mesh * 	mesh_,
						 const char * 		basename_,
						 cst_pI 			nrot_,
						 cst_pI 			axe_,
						 pI 			perm_,
						 pI                         rperm_,
						 cst_pI 			nvertex_axi_,
						 cst_pI                     nvertex_noaxi_);
  void 		print_scalar_extrude_axi	(const ns_mesh * 	mesh_,
						 const char * 		basename_,
						 cst_pI 			nrot_,
						 pI                       rperm_,
						 cst_pI 			nvertex_axi_,
						 cst_pI                   nvertex_noaxi_,
						 cst_pR                   scalar_);
  void 		print_vector_extrude_axi	(const ns_mesh * 	mesh_,
						 const char * 		basename_,
						 cst_pI 			nrot_,
						 cst_pI                   axe_,
						 pI                       rperm_,
						 cst_pI 			nvertex_axi_,
						 cst_pI                   nvertex_noaxi_,
						 cst_pR                   vectorx_,
						 cst_pR                   vectory_);

  void 		print_mesh_extrude_axi_permutation	(const ns_mesh * 	mesh_,
							 cst_pI 		axe_,
							 pI 			perm_,
							 pI 			rperm_,
							 pI                     nvertex_axi_,
							 pI                     nvertex_noaxi_);
  void 		print_mesh_extrude_axi_permutationP2	(const ns_mesh * 	mesh_,
							 cst_pI 		axe_,
							 pI 			perm_,
							 pI 			rperm_,
							 pI                   nvertex_axi_,
							 pI                   nvertex_noaxi_);

#if 0
  const nsCUBA_TRIA_INTERVAL_ST* 	ns_mesh_addcubature_interval	(ns_mesh * 		mesh_,
									 const I 		degree_);
  const nsCUBA_TRIA_INTERVAL_ST* 	ns_mesh_getcubature_interval	(const ns_mesh * 	mesh_,
									 const I 		degree_);
  const nsCUBA_TRIA_ST * 		ns_mesh_addcubature		(ns_mesh * 		mesh_,
									 const I 		degree_);
  const nsCUBA_TRIA_ST * 		ns_mesh_getcubature		(const ns_mesh * 	mesh_,
									 const I 		degree_);
#endif
  void ns_mesh_write_medit_with_codelm(const ns_mesh*s_,cst_pI codelm_,const char * name_,...);


#ifdef __cplusplus
}
#endif

#endif
