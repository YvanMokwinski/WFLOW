#ifndef __header_nsGLOBAL_h__
#define __header_nsGLOBAL_h__

#include "eStateSolver.h"
#include "eKindSystem.h"

#include "ns_constantes.h"
#include "ns_enum.h"
#include "eDim.h"
#include "ns_var.h"
#include "Cmdline.h"
#include "ens_method_transfert.h"

#include "Parameters.h"

#include "TransientScheme.h"
#include "Time.h"
#include "ns_config.h"
#include "ensBASIS.h"
#include "SparseStokes.h"
#include "Workelm.h"
/*#include "Variables.h"*/

#include "Fem.h"

#define NEW_SPARSE
#ifdef __cplusplus
extern "C"
{
#endif
  
  /* \brief structure for global variables */	
  typedef struct
  {
    int 			ithread;
    I 				iproc;
    void * 			mesh_usrptr;
    void * 			nonlinear_solver;
    void * 			transport_usrptr;

    /**
       \brief pointer to the finite element space of u 	
    */    
    pSpaceReadOnly 		space_u;
    /** 
	\brief pointer to the finite element space of p 	
    */
    pSpaceReadOnly 		space_p;
    /** 
	\brief pointer to the finite element space of dg	
    */
    pSpaceReadOnly 		space_dg;

    void * 			form_UP;
    
    I 				nelm;
    I 				nvertex;
    I 				nedge_boundary;
    I 				nedge;
    eMeshAdaptivity 		meshAdaptivityMethod;
    eTransientMethod 		method_transient_scheme_dg;


    pParameters 			parameters;

    Workelm			m_workelm;
    Time			m_timeInfo;
    TransientScheme		transient_schemes[__eKindEquation_ALL];

    CsfInfo 			m_csfInfo;
    Fem 			m_fem;

  
    pI				slip_iperm;
    /* SYSTEME NON LINEAIRE */ 
  
    /* SYSTEME LINEAIRE */ 
    SparseStokes*		sparseStokes;

    pSparse			sparseB;
    pSparse			sparseF;
    pSparse			sparseS;
    pSparse			sparseC;

    /*### 0/ METHODS ######################################## */
    cst_pR			jacelm;
    cst_pR			ddlcoo;
    cst_pI			ddlcod;
    cst_pI			ddlcnc;
    cst_pI			bpatch;
    cst_pI			patch ;
    cst_pI			adj ;
    cst_pR			jacedge;
    cst_pR			normaledge;
    cst_pR			trelm;
    const	char * 		filename;

    R				rv[__ens_rv_n];
    R 				ruser[16];

    /*### 1/ GENERAL VARIABLES ############################## */
    /* \brief error indicator 		*/
    /* \brief newton step counter		*/
    I 			newton_iter;
    /* \brief mesh element counter 	*/	
    I 			ielm;
    /* \brief adaptivity step counter 	*/	
    I 				adapt_iter;  
    /*### 2/ FREEDOM DEGREES ################################ */
    /* \brief total nddl */	
    I 				total_nddl;
    /* \brief total nddl for 1 comp. on u	*/	
    I 				nddlu;
    /* \brief nddl per elem u 		*/
    I 				nddlelmu;
    /* \brief total nddl 	for p 		*/	
    I 				nddlp;
    /* \brief nddl per elem for p		*/	
    I 				nddlelmp;
    /* \brief nddl per elm for dg		*/	
    I 				nddlelmdg;
    /* \brief total nddl for slip cond. 	*/	
    I 				nddl_slip;
    /* \brief shift to access to u	*/	
    I 				dec_ddlu;
    /* \brief shift to access to p	*/	
    I 				dec_ddlp;
    /* \brief shift to access to f	*/	

    I 				required_rwork_n;
    I 				required_iwork_n;
    I 				iwork_n;
    I 				blank_n;

    pI 				iwork;
    pI 				blank;
    ns_var 			var_u;
    ns_var 			var_p;
    char 				basename[512];
  } nsGLOBAL_ST,*RESTRICT pGlobal;

  typedef const nsGLOBAL_ST * RESTRICT pGlobalReadOnly;
  
  Err ns_build_residu	(pGlobal const self_,
			 pR 		global_rhs_,
			 cst_pR 	x_,
			 cst_pR 	xi_,
			 cst_pR 	xii_);

  Err ns_build_system	(pGlobal const self_,
			 cst_pI 	nx_,
			 cst_pR 	x_,
			 cst_pI 	xoff_);
  
  /* \brief structure for the Newton process */	
  typedef struct
  {
    Err (*build_residu)(pGlobal const self_,
			pR,
			cst_pR,
			cst_pR,
			cst_pR);
    Err (*build_system)(pGlobal const self_,
			cst_pI,
			cst_pR,
			cst_pI);
  } nsGLOBAL_NEWTON_DRIVER_ST;

  void	Global_initializePart1		(pGlobal 	const   self_,
					 ns_mesh *	const  	mesh_,
					 STR 			errmsg_,
					 Err*			err_);
  
  void Global_initializePart2		(pGlobal 	const   self_,
					 ns_mesh *	const  	mesh_,
					 const I 		n_,
					 const I 		ntimesteps_,
					 pR * 			x_,
					 STR 			errmsg_,
					 Err*			err_);



  pParametersReadOnly	GlobalReadOnly_get_Parameters	(pGlobalReadOnly 	const self_);

  pTimeReadOnly		GlobalReadOnly_get_Time		(pGlobalReadOnly 	const self_);
  pWorkelmReadOnly 	GlobalReadOnly_get_Workelm	(pGlobalReadOnly 	const self_);
  pFemReadOnly		GlobalReadOnly_get_Fem		(pGlobalReadOnly	const self_);

  pWorkelm 		Global_get_Workelm		(pGlobal 		const self_);

  pParameters 		Global_get_Parameters		(pGlobal 		const self_);

  pTime			Global_get_Time			(pGlobal 		const self_);

  pFem			Global_get_Fem			(pGlobal 		const self_);


  
  void 	ns_setelm			(pGlobal const 	self_,
					 const I 	jelm_,
					 cst_pR 	x_,
					 cst_pR 	xi_,
					 cst_pR 	xii_);
  
  
  Err sparse_A_init			(pGlobal const self_,
					 cst_eKindSystem kind_);
  
  Err sparse_F_init			(const I iproc_);

  Err sparse_B_init			(pGlobal const self_);
  
  void 	ns_build_system_divergence	(pGlobal const self_);
  
  void 	ns_build_system_slip		(pGlobal const self_,
					 cst_pI slip_iperm_,
					 cst_pI slip_nddlelm_);

  I 	ns_dg_fluxelm_boundary_get_coo	(pGlobal const  	G_,
					 const I 		jelm,
					 const I 		jadj,
					 pR 			pos_);
  void 	ns_dg_fluxelm_boundary		(pGlobal const  	G_,
					 cst_pR 		xa_,
					 const I 		jelm,
					 const I 		jadj,
					 pR 			rhs_,
					 cst_pR 		loc_,
					 cst_pR 		uvw_);


  void 			Global_main			(pGlobal 	const 	global_,
							 pCmdline 	const 	cmdline_,
							 STR 			errmsg_,
							 Err*			err_);

  void 			Global_precompute		(pGlobal	const 	self_,
							 ns_mesh * 	const 	mesh_,
							 STR 			errmsg_,
							 Err* 			err_);

  void 			Global_run			(pGlobal const 		self_,
							 eStateSolver 		job_[1],
							 cst_pI 		rwork_n_,
							 pR 			rwork_,
							 cst_pI 		iwork_n_,
							 pI 			iwork_);

  void 			Global_free			(pGlobal const  	global_);

  void 			nsGLOBAL_get_ddlcnc(pGlobal const G_,
					    const I _ielm,
					    pI _ddl);

  void 			nsGLOBAL_get_ddlcnc_u(pGlobal const G_,
					      const I _ielm,
					      pI _ddl);

  void 			nsGLOBAL_get_ddlcnc_v(pGlobal const 	G_,
					      const I 		_ielm,
					      pI 		_ddl);

  void 			nsGLOBAL_get_ddlcnc_p(pGlobal const 	G_,
					      const I 		_ielm,
					      pI 		_ddl);

  Err			Global_newton			(pGlobal const 		self_,
							 nsGLOBAL_NEWTON_DRIVER_ST * const driver_);


  L  			mkmake_ns_read_configfile	(const char * 	config_filename_,
							 R 		rinfo_[__ens_rinfo_ALL],
							 I  		iinfo_[__ens_iinfo_ALL],
							 STR 		sinfo_[__ens_sinfo_ALL],
							 L 		linfo_[__ens_linfo_ALL],
							 Err*		err_);

  Err 	nsGLOBAL_init_exact		(pGlobal const 	global_);
  Err 	nsGLOBAL_init_gauss		(pGlobal const 	global_,
					 const I 	degreelm_,
					 const I 	degreeface_,
					 cst_pI 	nf_,
					 cst_pI 	nu_,
					 cst_pI	np_);

  void 		nsGLOBAL_linsys_clear		(pGlobal const global_);

  /* \brief structure for triangle-based cubature		*/ 
  void 		nsGLOBAL_change_mesh		(pGlobal 	const  	self_,
						 const char * 	basename_);

  void 		nsGLOBAL_compute_volume		(pGlobal 	const self_);
  void 		nsGLOBAL_compute_cfl		(pGlobal 	const self_);
  void 		nsGLOBAL_transient_print	(pGlobalReadOnly 	const self_);

  Err 	nsGLOBAL_print_dgsol	(pGlobalReadOnly const self_,
				 const char * name_,
				 const ensBASIS shape_,
				 cst_pR x_,
				 cst_pI xoff_);
  Err 	nsGLOBAL_print_dgspace(pGlobalReadOnly const G_,
			       const char * 	name_,
			       const ensBASIS shape_);
  void 		nsGLOBAL_dg_print(pGlobalReadOnly const G_,
				  const char * 	name_,		 
				  const ensBASIS shape_,
				  cst_pR 	x_,
				  cst_pI 	xoff_);
  void 		nsGLOBAL_dg_print_sol(pGlobalReadOnly const G_,
				      const char * 	name_,		 	 
				      const ensBASIS shape_,
				      cst_pR 		x_,
				      cst_pI 		xoff_);


  void 		nsGLOBAL_dg_print_sol1	(pGlobalReadOnly const G_,
					 const char * 	name_,		 
					 const ensBASIS shape_,
					 cst_pR 		x_,
					 cst_pI 		xoff_);

  void 		nsGLOBAL_dg_print_mesh_interval(pGlobalReadOnly const G_,
						const char * 	name_,		 
						const ensBASIS shape_);


#if 0
  void 		LinsysVectors_update_previous_x(pGlobal 	self_,
						pLinsysVectors 	v_);
#endif  

#if 0
  pVariablesReadOnly	GlobalReadOnly_get_Variables	(pGlobalReadOnly	const self_);
  pVariables		Global_get_Variables		(pGlobal 		const self_);
#endif
  
#ifdef __cplusplus
}
#endif

#endif

