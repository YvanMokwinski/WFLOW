#ifndef __header_ns_var_h__
#define __header_ns_var_h__

#include "nsSPACE.h"
#include "ens_nrm.h"
#include "eTransientMethod.h"
#include "CsfInfo.h"
#include "Err.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /*\brief data structure for variables */
  typedef struct
  {
    ensBASIS 		lagr_id;
  
    L 			ortho;

    /* \brief the mesh 
     */  
    const ns_mesh * 	mesh;

    /* \brief the space 
     */  
    pSpaceReadOnly 	space;
  
    /* \brief The number of components 
     */
    I 			ncomp;

    /* \brief The number of dofs 
     */  
    I 			nddl;

  
    /*\brief number of freedom per elemenet */
    I 			nddlelm;
    /*\brief is discontinuous */
    I 			is_discont;
    /*\brief nddl if variable would be continuous, if discontinuous it is  */
    I         		cnddl;/*nb ddl si c'etait continu, = nddl dans le cas continu */
    /*\brief name of the variable */
    char 			name[1];
    /*\brief reference to the mesh */

    /*\brief x-component at the current time  */
    pR   			x;
    /*\brief x-component at the previous time  */
    pR   			xi;
    /*\brief x-component at the previous previous time  */
    pR   			xii;
  
    /*\brief y-component at the current time  */
    pR   			y;
    /*\brief y-component at the previous time  */
    pR   			yi;
    /*\brief y-component at the previous previous time  */
    pR   			yii;

    /*\brief z-component at the current time  */
    pR   			z;
    /*\brief z-component at the previous time  */
    pR   			zi;
    /*\brief z-component at the previous previous time  */
    pR   			zii;

    /*\brief adress of the current time */
    cst_pR 		t;
    /*\brief adress of the previous time */
    cst_pR 		ti;
    /*\brief adress of the previous previous time */
    cst_pR 		tii;
  
    R 			nrms	[__emk_nrm_n];
    R 			nrms_i	[__emk_nrm_n];
    R 			nrms_ii	[__emk_nrm_n];

    R 			analytic_nrms	[__emk_nrm_n];
    R 			analytic_nrms_i	[__emk_nrm_n];
    R 			analytic_nrms_ii[__emk_nrm_n];

  } ns_var;

  /* ############################################################################ */

  void ns_var_cpy_xi2xii(ns_var * self_);

#ifndef NDEBUG

  const ns_mesh * 	ns_var_mesh		(const ns_var*const F_);
  pSpaceReadOnly 		ns_var_space		(const ns_var*const F_);

  ensBASIS		ns_var_shape		(const ns_var*const F_);
  L			ns_var_isdiscont	(const ns_var*const F_);
  L			ns_var_isortho		(const ns_var*const F_);
  I     			ns_var_nddlelm		(const ns_var*const F_);
  I     			ns_var_ncomp		(const ns_var*const F_);
  I     			ns_var_nddl		(const ns_var*const F_);

  cst_pR 			ns_var_t		(const ns_var*const F_);
  cst_pR 			ns_var_ti		(const ns_var*const F_);
  cst_pR 			ns_var_tii		(const ns_var*const F_);

  cst_pR 			cst_ns_var_x		(const ns_var*const F_);
  cst_pR 			cst_ns_var_xi		(const ns_var*const F_);
  cst_pR 			cst_ns_var_xii		(const ns_var*const F_);

  cst_pR 			cst_ns_var_y		(const ns_var*const F_);
  cst_pR 			cst_ns_var_yi		(const ns_var*const F_);
  cst_pR 			cst_ns_var_yii		(const ns_var*const F_);

  cst_pR 			cst_ns_var_z		(const ns_var*const F_);
  cst_pR 			cst_ns_var_zi		(const ns_var*const F_);
  cst_pR 			cst_ns_var_zii		(const ns_var*const F_);

  pR 			ns_var_x		(ns_var*const F_);
  pR 			ns_var_xi		(ns_var*const F_);
  pR 			ns_var_xii		(ns_var*const F_);

  pR 			ns_var_y		(ns_var*const F_);
  pR 			ns_var_yi		(ns_var*const F_);
  pR 			ns_var_yii		(ns_var*const F_);

  pR 			ns_var_z		(ns_var*const F_);
  pR 			ns_var_zi		(ns_var*const F_);
  pR 			ns_var_zii		(ns_var*const F_);


#else

#define         ns_var_isdiscont(_v)	(((_v)->is_discont)?__emk_yes:__emk_no)
#define         ns_var_isortho(_v)	(_v)->ortho
#define         ns_var_t(_v)		(_v)->t
#define         ns_var_ti(_v)		(_v)->ti
#define         ns_var_tii(_v)		(_v)->tii

#define         ns_var_x(_v)		(_v)->x
#define         ns_var_xi(_v)		(_v)->xi
#define         ns_var_xii(_v)		(_v)->xii

#define         ns_var_y(_v)		(_v)->y
#define         ns_var_yi(_v)		(_v)->yi
#define         ns_var_yii(_v)		(_v)->yii

#define         ns_var_z(_v)		(_v)->z
#define         ns_var_zi(_v)		(_v)->zi
#define         ns_var_zii(_v)		(_v)->zii

#define         cst_ns_var_x(_v)	(_v)->x
#define         cst_ns_var_xi(_v)	(_v)->xi
#define         cst_ns_var_xii(_v)	(_v)->xii

#define         cst_ns_var_y(_v)	(_v)->y
#define         cst_ns_var_yi(_v)	(_v)->yi
#define         cst_ns_var_yii(_v)	(_v)->yii

#define         cst_ns_var_z(_v)	(_v)->z
#define         cst_ns_var_zi(_v)	(_v)->zi
#define         cst_ns_var_zii(_v)	(_v)->zii

#define         ns_var_shape(_v)	(_v)->lagr_id
#define         ns_var_space(_v)	(_v)->space
#define         ns_var_mesh(_v)		(_v)->mesh
#define 	ns_var_nddlelm(_v) 	(_v)->nddlelm
#define 	ns_var_nddl(_v) 	(_v)->nddl
#define         ns_var_ncomp(_v)	(_v)->ncomp

#endif



  ns_var * 	ns_mesh_killvar			(ns_mesh*		mesh_,
						 ns_var*		v_);

  ns_var * 	ns_mesh_addvector		(ns_mesh * 		mesh_,
						 const ensBASIS  	shape_,
						 const L 		ortho_);

  ns_var * 	ns_mesh_addscalar		(ns_mesh * 		mesh_,
						 const ensBASIS  	shape_,
						 const L 		ortho_);
  ns_var * 	ns_mesh_addtensor		(ns_mesh * 		mesh_,
						 const ensBASIS  	shape_,
						 const L 		ortho_);

  /* ############################################################################
     transfert entre deux maillages */
  void* 		ns_var_transfert		(ns_var * 		self_,
							 const ns_var * 	oldVariable_,
							 const I 		method_,
							 void * 		old_localization_,
							 Err*			err_);


  void 		ns_var_free			(ns_var*		self_);
  ns_var*		ns_var_kill			(ns_var*		self_);
  void 		ns_var_clear			(ns_var*		self_);
  void 		ns_var_info			(const ns_var * 	self_);


  /***** UTILS BASIC */

  void 		ns_var_zero			(ns_var * 			self_);
  void 		ns_var_zeroi			(ns_var * 			self_);
  void 		ns_var_zeroii			(ns_var * 			self_);
  void 		ns_var_lagr2l2ortho		(ns_var * 			self_);
  void 		ns_var_lagr2l2ortho_i		(ns_var * 			self_);
  void 		ns_var_lagr2l2ortho_ii		(ns_var * 			self_);
  void 		ns_var_l2ortho2lagr		(ns_var *			self_);


  /***** UTILS I/O */
  void 		ns_var_write_medit		(const ns_var * 		self_,
						 const char * 			name_,
						 ...);

  Err 	ns_var_read			(ns_var * 	self_,
					 const char* 	name_,
					 const I 	i_);

  Err 	ns_var_read_medit		(ns_var * 	self_,
					 const char* 	name_);
  void 		ns_var_write_medit_with_mesh	(const ns_var * f_,
						 const char * 	name_,
						 ...);
  Err 	ns_var_read_bin			(ns_var * 	self_,
					 const char* 	name_);
  Err 	ns_var_write_bin		(const ns_var * self_,
					 const char* 	name_);

  void 		ns_var_extrude_time_medit	(const ns_var * self_,
						 const char * 	name_,
						 const I 	nframe_,
						 const I 	nz_,
						 cst_pR 	z_,
						 const R 	coeff_,
						 cst_pR        *scalars_);

  void 		ns_var_print			(const ns_var * 	self_,
						 const char * 		name_,
						 cst_pI 		itime_interval_,
						 cst_pI 		itime_);



  /***** UTILS CSF */
  void 		ns_var_csf			(const ns_var * 	f_,
						 ns_var * 		nabla_,
						 ns_var*		tension_,
						 ns_var * 		normal_,
						 ns_var*		kappa_,
						 pCsfInfoReadOnly 	csf_);

  /***** UTILS FEM */
  void 		ns_var_expand			(const ns_var* c_,ns_var * d_);
  void 		ns_var_reduce			(const ns_var* d_,ns_var * c_);
  void 		ns_var_reduce_nosum		(const ns_var* d_,ns_var * c_);


  /***** UTILS NRMS */
  void 		ns_var_nrmjump			(const ns_var * self_,
						 R nrms_[6]);

  void 		ns_var_nrm_all			(const ns_var*	self_,
						 pR 		nrm_);

  /***** UTILS DERIVATIVE */
  void  		ns_var_nabla_scalar		(ns_var*			nabla_,
							 const ns_var * 		f_,
							 cst_nsOP_INFO			op_info_);

  void 		ns_var_kappa			(const ns_var * 		self_,
						 ns_var*			divergence_,
						 cst_nsOP_INFO			op_info_,
						 const L 			filtering_,
						 cst_pR 			tol_);

  void 		ns_var_normalize		(ns_var * 			self_,
						 cst_nsOP_INFO 			op_info_,
						 cst_pR 			tol_);

  void  		ns_var_scalar_lsq		(ns_var * 			f_,
							 cst_eSmoothingMethod		method_,
							 const I 			increment_degree_);



  void 		ns_var_smoothl2ortho		(const ns_var*			self_,
						 ns_var*			w_,
						 const ensBASIS 		basis_shape_);





  void 		ns_var_cont2discont		(const ns_var*			self_,
						 ns_var*			w_);






  void  		ns_var_divergence_vector		(ns_var*			divergence_,
								 const ns_var * 		self_,
								 cst_nsOP_INFO			op_info_,
								 cst_pI			selection_n_,
								 cst_pI			selection_);

  void 		ns_var_vector_scalarxvector_dirac	(ns_var * 			self_,
							 cst_eOperator 		op_,
							 const ns_var* 			s_,
							 const ns_var* 			n_,
							 const ns_var* 			dist_,
							 pSmoothedDiracReadOnly		smoothedDiracReadOnly_);




  ns_var * 	ns_var_new_malloc		(const ns_mesh * 		mesh_,
						 const char * 			name_,
						 const ensBASIS 		shape_,
						 const L           	is_discont_,
						 const I            	ncomp_,
						 cst_pR 			t_,
						 cst_pR 			ti_,
						 cst_pR 			tii_);




  void 		ns_var_def			(ns_var * 		self_,
						 const ns_mesh * 	mesh_,
						 const char * 		name_,
						 const ensBASIS 	shape_,
						 const L           is_discont_,
						 cst_pI                ncomp_,
						 cst_pR 		t_,
						 cst_pR 		ti_,
						 cst_pR 		tii_,
						 pR 			x_,
						 pR 			xi_,
						 pR 			xii_);

  void 		ns_var_def_malloc		(ns_var * 		self_,
						 const ns_mesh * 	mesh_,
						 const char * 		name_,
						 const ensBASIS 	shape_,
						 const L           is_discont_,
						 cst_pI                ncomp_,
						 cst_pR 		t_,
						 cst_pR 		ti_,
						 cst_pR 		tii_);

  void 		ns_var_d2c			(const ns_var * self_,
						 ns_var *  	f_);

  void 		ns_var_cpy2xi			(ns_var * 	self_,/* SI LA VARIABLEEST CONTINU */
						 cst_pR 	a_,/*a continu*/
						 cst_pI 	aoff_);
  void 		ns_var_copyx			(ns_var * 	self_,
						 const ns_var * v2_);

  void 		ns_var_cpy_xi2x			(ns_var * 	self_);
  void 		ns_var_cpy_x2xi			(ns_var * 	self_);

  void 		ns_var_welm			(const ns_var*	self_,
						 cst_pI 	ielm_,
						 pR 		w_,
						 cst_pI 	woff_);

  void 		ns_var_welmi			(const ns_var*	self_,
						 cst_pI 	ielm_,
						 pR 		wi_,
						 cst_pI 	woff_);

  void 		ns_var_welmii			(const ns_var*	self_,
						 cst_pI 	ielm_,
						 pR 		wii_,
						 cst_pI 	woff_);

  void 		ns_var_welms			(const ns_var*	self_,
						 cst_pI 	nelm_,
						 cst_pI 	elm_,
						 pR 		w_,
						 cst_pI 	woff_);

  void 		ns_var_nabla_shape		(const ns_var*	self_,
						 ns_var*	nabla_);

  void		ns_var_divergence_shape		(const ns_var*	self_,
						 ns_var*	divergence_);

  ns_var*		ns_var_nabla_shape_new		(const ns_var*v_);


  /*-----------------------------------------------------------------------------*/

  /*-----------------------------------------------------------------------------*/



  /*-----------------------------------------------------------------------------*/

  /*-----------------------------------------------------------------------------*/

  /*-----------------------------------------------------------------------------*/

  /*-----------------------------------------------------------------------------*/

#if 0

  void 		ns_var_csf_color_shape_ph	(ns_var*			csf_,
						 const ns_var*			v_,
						 ns_var*			nabla_,
						 ns_var*			normal_,
						 ns_var*			kappa_,
						 const nsCUBA_TRIA_ST*	cubature_);

  void 		ns_var_normalization_ph	(ns_var * 		normal_,
					 const ns_var *		v_,
					 const nsCUBA_TRIA_ST*	cubature_,
					 cst_pR 		tol_);


  void 		ns_var_vector_scalarxvector_ph	(ns_var * 		vector_,
						 const ns_var * 	scalar_,
						 const ns_var * 	v_,
						 const nsCUBA_TRIA_ST*	cubature_);

  void  		ns_var_nabla_vector_lsq	(const ns_var * 		vec_,
						 const nsCUBA_TRIA_ST * 	cubature_,
						 cst_eSmoothingMethod		method_,
						 const I 			increment_degree_,
						 ns_var*			v_);
  void  		ns_var_nabla_scalar_lsq	(const ns_var * 		f_,				
						 cst_eSmoothingMethod		method_,
						 const I 			increment_degree_,
						 ns_var*			v_);
  void  		ns_var_f_lsq		(const ns_var * 		f_,
						 const nsCUBA_TRIA_ST *	cubature_,
						 cst_ensHEAVISIDE		heaviside_,
						 cst_ensDIRAC			dirac_,
						 cst_pR 			eps_,
						 cst_eSmoothingMethod 		method_,
						 const I			increment_degree_,
						 ns_var*			v_);
  mkZCURVE ns_var_iso		(const ns_var *                	self_,
				 cst_nsISO_INFO			iso_info_);

  void 	ns_var_iso_ih 		(ns_var * 			self_,
				 mkZCURVE  			iso_);
  void 	ns_var_iso_ph 		(ns_var * 			self_,
				 mkZCURVE  			iso_);
#endif

#ifdef __cplusplus
}
#endif

#endif
