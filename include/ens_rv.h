#ifndef __header__ens_rv_h__
#define __header__ens_rv_h__


enum __ens_rv 	{__ens_rv_iweber=0,
		 __ens_rv_ireynold,
		 __ens_rv_ifroude,
		 __ens_rv_vmin,
		 __ens_rv_vmax,
		 __ens_rv_coeff_ratio_viscosity,
		 __ens_rv_coeff_ratio_density,
		 __ens_rv_dt,
		 __ens_rv_dti,
		 __ens_rv_idt,
		 __ens_rv_midt,
		 __ens_rv_ratio_dti,
		 __ens_rv_iratio_dti,
		 __ens_rv_eps,
		 __ens_rv_ieps,
		 __ens_rv_n };

typedef enum __ens_rv ens_rv;

typedef const ens_rv cst_ens_rv;

#endif
