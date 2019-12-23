#ifndef __header__ens_rres_h__
#define __header__ens_rres_h__

enum __ens_rres { __ens_rres_error=0,
		  __ens_rres_t,
		  __ens_rres_dt,
		  __ens_rres_nelm,
		  __ens_rres_nvertex,
		  __ens_rres_nddlu,
		  __ens_rres_nddlp,
		  __ens_rres_nddlf,
		  __ens_rres_uL1,
		  __ens_rres_uL2,
		  __ens_rres_uLInf,
		  __ens_rres_uH1,
		  __ens_rres_pL1,
		  __ens_rres_pL2,
		  __ens_rres_pLInf,
		  __ens_rres_pH1,
		  __ens_rres_fL1,
		  __ens_rres_fL2,
		  __ens_rres_fLInf,
		  __ens_rres_fH1,
		  __ens_rres_area,
		  __ens_rres_length_freesurface,		  
		  __ens_rres_etime_total,
		  __ens_rres_etime_loading,
		  __ens_rres_etime_init_data,
		  __ens_rres_etime_kill_data,
		  __ens_rres_etime_algo_diff,
		  __ens_rres_etime_algo_normal,
		  __ens_rres_etime_algo_kappa,
		  __ens_rres_etime_algo_linsol,		       		       		      		       
		  __ens_rres_n };

typedef enum __ens_rres ens_rres;
typedef const ens_rres cst_ens_rres;

#endif
