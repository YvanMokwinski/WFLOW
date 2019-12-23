#ifndef __header_DgTransport_h__
#define __header_DgTransport_h__
#include <math.h>
#include "Blas.h"
#include "Type.h"
#include "mkS.h"

#define DG_CONSERVATIVE_WEAK_FORM 	NO

#define DGERR_MEMORY  2
#define DGERR_USER    3

typedef enum __eDgInfo
  {

    DG_IA_lc=0,
    DG_IA_lc_elm=1,
    DG_IA_lc_face=2,
    DG_IA_lc_elm_a=3,
    DG_IA_lc_u=4,
    DG_IA_lc_elm_u0=5,
    DG_IA_lc_elm_u1=6,
    DG_I_nTot=11,
    DG_I_nTotElm=12,
    DG_I_nTotFace=13,
    DG_RA_lc=14,
    DG_I_lc_len=16,
    DG_RA_bmat=17,
    DG_RA_bmatx=18,
    DG_I_bmat_n=19,
    DG_I_bmat_m=20,
    DG_I_bmat_len=21,
    DG_ERR=22,
    DG_RA_EVAL_TETA_U=23,

    DG_I_QFACE_N=28,
    DG_I_QELM_N=29,

    DG_I_TRIAL_FAMILY=30,
    DG_I_TRIAL_DEGREE=31,
    DG_I_TRIAL_NBASIS=32,

    
    DG_I_TEST_FAMILY=33,
    DG_I_TEST_DEGREE=34,
    DG_I_TEST_NBASIS=35,
    
    DG_I_TETA_FAMILY=36,
    DG_I_TETA_DEGREE=37,
    DG_I_TETA_NBASIS=38,
    
    DG_I_TETA_U_FAMILY=39,
    DG_I_TETA_U_DEGREE=40,
    DG_I_TETA_U_NBASIS=41,
    
    DG_I_TETA_A_FAMILY=42,
    DG_I_TETA_A_DEGREE=43,
    DG_I_TETA_A_NBASIS=44,

    DG_I_n=64

  } eDgInfo;


#define DG_r_lc        			0
#define DG_r_n         			1

#define DG_ires_err 			0 
#define DG_ires_convergence 		1
#define DG_ires_iter_gauss_seidel 	2
#define DG_ires_required_iw_n 		3
#define DG_ires_required_rw_n 		4
#define DG_ires_n 			5

#define DG_rres_max 			0 
#define DG_rres_nrmL2 			1 
#define DG_rres_nrmLInf			2 
#define DG_rres_areaL1 			3 
#define DG_rres_jumpL2 			4 
#define DG_rres_johnson			5 
#define DG_rres_n 			6

#ifdef __cplusplus
extern "C"
{
#endif

  //
  //
  //
  void  dgadvection	(cst_mkS 	s_teta_a_,
			 cst_mkS 	s_teta_u_,
			 cst_mkS 	s_teta_,
			 cst_mkS 	s_test_,
			 cst_mkS 	s_trial_,
			 cst_pI         iinfo_n_,
			 pI		iinfo_,
			 cst_pI		rinfo_n_,
			 pR 		rinfo_,
			 cst_pI 	rwork_n_,
			 pR		rwork_);
  
  void DgTransport_compute(cst_pR	xu_,
			   cst_pR	rhs_,
			   cst_pI	rhsoff_,
			   cst_pR	data_u_,
			   cst_pI	data_uoff_,
			   cst_pR	data_v_,
			   cst_pI	data_voff_,
			   cst_pR	sol_,
			   cst_pI	soloff_,
			   pR 		corr_,
			   cst_pI 	corroff_ ,
			   cst_pR 	t_,
			   
			   cst_pI 	nvertex_,
			   cst_pI 	nelm_,
			   cst_pR 	coo_,
			   cst_pI 	cooff_,
			   cst_pI 	cnc_,
			   cst_pI 	cncoff_,
			   cst_pI 	adj_,
			   cst_pI 	adjoff_,
			   cst_pI 	vcod_,
			   cst_pI 	noboundary_cod_,


			 
			   // cst_mkZ   	zone_,
			   cst_pI 	rwork_n_,
			   pR  		rwork_,
			   cst_pI 	iwork_n_,
			   pI		iwork_,
			   cst_pR 	rinfo_,
			   const I  	iinfo_[DG_I_n],
			   R 	rres_[DG_rres_n],
			   I 		ires_[DG_ires_n]);

  void dgadvection_solve(cst_pR 	xa_,
			 cst_pR 	xu_,
			 cst_pR		rhs_,
			 cst_pI  	rhsoff_,
			 cst_pI         cnc_u_,
			 cst_pI         cncoff_u_,
			 cst_pR 	data_u_,
			 cst_pR		data_v_,
			 cst_pR		sol_,
			 cst_pI 	soloff_,
			 pR		corr_,
			 cst_pI 	corroff_ ,
			 cst_pR 	t_,
			 cst_pI 	nelm_,
			 cst_pR 	coo_,
			 cst_pI 	cooff_,
			 cst_pI 	cnc_,
			 cst_pI 	cncoff_,
			 cst_pI 	adj_,
			 cst_pI 	adjoff_,
			 cst_pI 	vcod_,
			 cst_pI 	noboundary_cod_,
			 
			 cst_pI 	rwork_n_,
			 pR  		rwork_,
			 cst_pI 	iwork_n_,
			 pI		iwork_,
			 cst_pR 	rinfo_,
			 const I 	iinfo_[DG_I_n],
			 R 	 	rres_[DG_rres_n],
			 I 	 	ires_[DG_ires_n]);
  
  void DgTransport_routine(cst_pR 	xa_,
			   cst_pR 	xu_,
			   cst_pR		rhs_,
			   cst_pI  	rhsoff_,
			   cst_pI         cnc_u_,
			   cst_pI         cncoff_u_,
			   cst_pR 	data_u_,
			   cst_pR		data_v_,
			   cst_pR		sol_,
			   cst_pI 	soloff_,
			   pR		corr_,
			   cst_pI 	corroff_ ,
			
			   cst_pR 	t_,
			
			   cst_pI 	nelm_,
			   cst_pR 	coo_,
			   cst_pI 	cooff_,
			   cst_pI 	cnc_,
			   cst_pI 	cncoff_,
			   cst_pI 	adj_,
			   cst_pI 	adjoff_,
			   cst_pI 	vcod_,
			 cst_pI 	noboundary_cod_,
			
			 cst_pI 	rwork_n_,
			 pR  		rwork_,
			 cst_pI 	iwork_n_,
			 pI		iwork_,
			 cst_pR 	rinfo_,
			 const I 	iinfo_[DG_I_n],
			 R 	 	rres_[DG_rres_n],
			 I 	 	ires_[DG_ires_n]);


#if 0
void  DgTransport_init(cst_mkS 		s_teta_a_,
		       cst_mkS 		s_teta_u_,
		       cst_mkS 		s_teta_,
		       cst_mkS 		s_test_,
		       cst_mkS 		s_trial_,
		       I		iinfo_[DG_I_n],
		       cst_pI		rinfo_n_,
		       pR 		rinfo_,
		       cst_pI 		rwork_n_,
		       pR		rwork_);
#endif



typedef struct
{
  I 		dg_iinfo[DG_I_n];
  R		dg_rres	[DG_rres_n];
  I			dg_ires	[DG_ires_n];
  I 		dg_rinfo_n;
  pR 			dg_rinfo;
  I 		dg_rwork_n;
  pR 			dg_rwork;

  I 		dgerr_iinfo[DG_I_n];
  R		dgerr_rres	[DG_rres_n];
  I			dgerr_ires	[DG_ires_n];
  I 		dgerr_rinfo_n;
  pR 			dgerr_rinfo;
  I 		dgerr_rwork_n;
  pR 			dgerr_rwork;

  I 		iwork_n;
  pI 			iwork;

  pR 			f;
  pR 			fi;
  pR 			fii;
  pR 			fiii;
  pR 			rhs;
  pR 			corr;

  pR 			e;
  pR 			ei;
  pR 			eii;
  pR 			eiii;
  pR 			erhs;
  pR 			ecorr;

} mkDG_COMMON_ST,*mkDG_COMMON;

typedef const mkDG_COMMON_ST cst_mkDG_COMMON_ST;
typedef cst_mkDG_COMMON_ST*cst_mkDG_COMMON;

size_t mkDG_COMMON_fwrite	(cst_mkDG_COMMON 	dg_common_,
				 FILE*			fich);

size_t mkDG_COMMON_fread	(mkDG_COMMON 	 	dg_common_,
				 FILE*			fich);

#ifdef __cplusplus
}
#endif


#endif
