#ifndef __header_nsGLOBAL_h__
#define __header_nsGLOBAL_h__

#include "eStateSolver.h"
#include "eKindSystem.h"
#include "ns_constantes.h"
#include "ns_enum.h"
#include "eDim.h"

#include "ns_var.h"

#if 0
#include "nsMETRIC.h"
#endif
#include "Cmdline.h"
#include "ens_method_transfert.h"
#include "Parameters.h"
#include "TransientScheme.h"
#include "Time.h"
#include "ns_config.h"

#include "ensBASIS.h"
#include "Linsys.h"
#include "LinsysVectors.h"
#include "SparseStokes.h"
#define NEW_SPARSE
Err ns_build_residu	(const I 	iproc_,
			 pR 		global_rhs_,
			 cst_pR 	x_,
			 cst_pR 	xi_,
			 cst_pR 	xii_);
Err ns_build_system	(const I 	iproc_,
			 cst_pI 	nx_,
			 cst_pR 	x_,
			 cst_pI 	xoff_);

/* \brief structure for the Newton process */	
typedef struct
{
  Err (*build_residu)(const I,pR ,cst_pR ,cst_pR ,cst_pR);
  Err (*build_system)(const I,cst_pI ,cst_pR ,cst_pI );
} nsGLOBAL_NEWTON_DRIVER_ST;

typedef struct
{
  pR 		p;
  pR 		pi;
  pR 		pii;
  pR 		u;
  pR 		ui;
  pR 		uii;
  pR 		v;
  pR 		vi;
  pR 		vii;
  pR 		w;
  pR 		wi;
  pR 		wii;
} nsHANDLE_VAR_ST;



Err nsHANDLE_VAR_def(nsHANDLE_VAR_ST *const hdl_,
		     const I iproc_,
		     pLinsysVectors linsys_vectors_,
		     cst_pI nddlu_,
		     cst_pI dec_ddlu_,			 
		     cst_pI dec_ddlp_);

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
}  nsWORKELM_ST;





#define perf_t0 		0
#define perf_t_residu 		1
#define perf_t_sys 		2
#define perf_t_solve 		3
#define perf_tot_residu 	4
#define perf_tot_sys 		5
#define perf_tot_solve 		6
#define perf_t_work		7
#define perf_t_iter 		8
#define perf_tot_iter 		9

typedef struct
{
  R 	nrmr	[32];
  R 	nrmc	[32];
  R 	nrmh	[32];
  R 	nrmcs	[32];
  R 	nrmcu	[32];
  R 	nrmcp	[32];
  R 	nrmrs	[32];
  R 	nrmru	[32];
  R 	nrmrp	[32];
} nsGLOBAL_NEWTON_RESULTS_ST;

void nsGLOBAL_NEWTON_RESULTS_def(nsGLOBAL_NEWTON_RESULTS_ST*const v_);

typedef struct
{
  R p[32];
} nsGLOBAL_NEWTON_PERF_ST;

void nsGLOBAL_NEWTON_PERF_def	(nsGLOBAL_NEWTON_PERF_ST*const perf_);

#define _build_partialoff 		( _nu*_nu*_nu*_ndg * 3 + _nu*_ndg*_ndg*4   +  _nu*_np   )
#define _build_twice_partialoff  	( _nu*_nu*_ndg*2+_nu*_ndg*_ndg*2)

typedef struct
{
  eElement elm;
 
  /* L ORDRE EST IMPORTANT */
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
  
} nsFEM_DATA_ST;


void 	nsFEM_DATA_setelm		(nsFEM_DATA_ST*const 	fem_,
					 nsWORKELM_ST*const 	elm_);
Err nsFEM_DATA_def			(nsFEM_DATA_ST*const 	fem_,
					 cst_eElement 		elm_);

Err nsFEM_DATA_def_from_quadrature	(nsFEM_DATA_ST*const 	fem_,
					 cst_pI		qn_,
					 cst_pR		qw_,
					 cst_pR		qp_,
					 cst_pI		qoff_);


/* \brief structure for global variables */	
typedef struct
{
  int 				ithread;
  I 			iproc;
  void * 			mesh_usrptr;
  void * 			nonlinear_solver;
  void * 			transport_usrptr;
  /* \brief pointer to the finite element space of u 	*/
  const nsSPACE_ST * 		space_u;
  /* \brief pointer to the finite element space of p 	*/
  const nsSPACE_ST * 		space_p;
  /* \brief pointer to the finite element space of dg	*/
  const nsSPACE_ST * 		space_dg;
  void * form_UP;
  I 				nelm;
  I 				nvertex;
  I 				nedge_boundary;
  I 				nedge;
  eMeshAdaptivity 		meshAdaptivityMethod;
  eTransientMethod 		method_transient_scheme_dg;

  Parameters 			parameters;
  nsWORKELM_ST			workelm;
  Time 				timeInfo;
  TransientScheme		transient_schemes[__eKindEquation_ALL];
  nsHANDLE_VAR_ST 		hdl;
  CsfInfo 			csfInfo;

  pI				slip_iperm;
  /* SYSTEME NON LINEAIRE */ 
  nsGLOBAL_NEWTON_PERF_ST 	newton_perf;
  nsGLOBAL_NEWTON_RESULTS_ST 	newton_results;
  /* SYSTEME LINEAIRE */ 



  pLinsysVectors		linsysVectors;
  SparseStokes			sparseStokes;
  pSparse			sparseB;
  pSparse			sparseF;
  pSparse			sparseS;
  pSparse			sparseC;
  pLinsys			linsysA;
  pLinsys			linsysF;
  pLinsys			linsysB;
  pLinsys			linsysS;
  pLinsys			linsysC;
  nsFEM_DATA_ST 		fem;
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
  const	char * 			filename;
  R			rv[__ens_rv_n];
  R 			ruser[16];

  /*### 1/ GENERAL VARIABLES ############################## */
  /* \brief error indicator 		*/
  /* \brief newton step counter		*/
  I 			newton_iter;
  /* \brief mesh element counter 	*/	
  I 			ielm;
  /* \brief adaptivity step counter 	*/	
  I 			adapt_iter;  
  /*### 2/ FREEDOM DEGREES ################################ */
  /* \brief total nddl */	
  I 			total_nddl;
  /* \brief total nddl for 1 comp. on u	*/	
  I 			nddlu;
  /* \brief nddl per elem u 		*/
  I 			nddlelmu;
  /* \brief total nddl 	for p 		*/	
  I 			nddlp;
  /* \brief nddl per elem for p		*/	
  I 			nddlelmp;
  /* \brief nddl per elm for dg		*/	
  I 			nddlelmdg;
  /* \brief total nddl for slip cond. 	*/	
  I 			nddl_slip;
  /* \brief shift to access to u	*/	
  I 			dec_ddlu;
  /* \brief shift to access to p	*/	
  I 			dec_ddlp;
  /* \brief shift to access to f	*/	

  I 			required_rwork_n;
  I 			required_iwork_n;
  I 			iwork_n;
  I 			blank_n;

  pI 				iwork;
  pI 				blank;
  ns_var 			var_u;
  ns_var 			var_p;
#if 0
  nsMETRIC_ST* 			metric;
  nsMETRIC_ST* 			metric_end_time_interval;
#endif
  char 				basename[512];
} nsGLOBAL_ST,Global,*RESTRICT pGlobal;

typedef const Global * RESTRICT cst_pGlobal;

pParameters 	Global_get_Parameters		(pGlobal self_);
cst_pParameters	cst_Global_get_Parameters	(cst_pGlobal self_);

pTime		Global_get_Time		(pGlobal self_);
cst_pTime	cst_Global_get_Time	(cst_pGlobal self_);


void 	ns_setelm			(const I iproc_,
					 const I jelm_,
					 cst_pR x_,
					 cst_pR xi_,
					 cst_pR xii_);


Err sparse_A_init			(const I 	iproc_,
					 cst_eKindSystem kind_);
Err sparse_F_init			(const I iproc_);
Err sparse_B_init			(const I iproc_);
void 	ns_build_system_divergence	(const I iproc_);
void 	ns_build_system_slip		(const I iproc_,
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





Err 	   	nsGLOBAL_initialize		(pGlobal 	const   G_,
							 ns_mesh *	const  	mesh_);

Err 		Global_main			(pGlobal 	const 	global_,
							 pCmdline 	const 	cmdline_);

Err 		Global_precompute		(pGlobal const 	global_,
							 ns_mesh*const 		mesh_);

Err 		nsGLOBAL_init			(pGlobal const  	global_,
							 const char * 		filename_,
							 const char * 		basename_);


void 			Global_run			(pGlobal const 		self_,
							 eStateSolver 	job_[1],
							 cst_pI 		rwork_n_,
							 pR 			rwork_,
							 cst_pI 		iwork_n_,
							 pI 			iwork_);

void 			Global_free			(pGlobal const  	global_);

void nsGLOBAL_get_ddlcnc(pGlobal const G_,
			 const I _ielm,
			 pI _ddl);
void nsGLOBAL_get_ddlcnc_u(pGlobal const G_,
			 const I _ielm,
			 pI _ddl);
void nsGLOBAL_get_ddlcnc_v(pGlobal const G_,
			 const I _ielm,
			 pI _ddl);
void nsGLOBAL_get_ddlcnc_p(pGlobal const G_,
			 const I _ielm,
			 pI _ddl);

Err			Global_newton			(pGlobal const 		self_,
							 nsGLOBAL_NEWTON_DRIVER_ST * const driver_);


L  	mkmake_ns_read_configfile	(const char * 		config_filename_,
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
void 		nsGLOBAL_transient_print	(cst_pGlobal 	const self_);

Err 	nsGLOBAL_print_dgsol	(cst_pGlobal const self_,
					 const char * name_,
					 const ensBASIS shape_,
					 cst_pR x_,
					 cst_pI xoff_);
Err 	nsGLOBAL_print_dgspace(cst_pGlobal const G_,
				       const char * 	name_,
				       const ensBASIS shape_);
void 		nsGLOBAL_dg_print(cst_pGlobal const G_,
				  const char * 	name_,		 
				  const ensBASIS shape_,
				  cst_pR 	x_,
				  cst_pI 	xoff_);
void 		nsGLOBAL_dg_print_sol(cst_pGlobal const G_,
				      const char * 	name_,		 	 
				      const ensBASIS shape_,
				      cst_pR 		x_,
				      cst_pI 		xoff_);


void 		nsGLOBAL_dg_print_sol1	(cst_pGlobal const G_,
					 const char * 	name_,		 
					 const ensBASIS shape_,
					 cst_pR 		x_,
					 cst_pI 		xoff_);
void 		nsGLOBAL_dg_print_mesh_interval(cst_pGlobal const G_,
						const char * 	name_,		 
						const ensBASIS shape_);


#endif

