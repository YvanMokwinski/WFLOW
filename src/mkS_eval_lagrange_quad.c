#include "mkSYS.h"
#include "mkS.h"
#include "mok_m.h"
#include "mkS_eval.h"


#if __mk_debug__
void mkS_lagrange_localspl_quad( cst_mpsint 	degree,
				 mkR		coo_,
				 cst_mkI	off_)
{
  __mmkS_lagrange_localspl_quad(degree,coo_,off_);
}
#endif


void mkS_dy_lagrange_quad(	     cst_mpsint degree,
					     cst_mpsint n,
					     mpsreal r,
					     cst_mpsint roff_,
					     cst_mpsreal  p,
					     cst_mpsint poff_,
					     mpsreal rwork,
					     cst_mpsint rwork_n,
					     mpsint err_)
{
  nsINT info,nshapes = (degree[0]+1)*(degree[0]+1);
  nsINT perm[128];
  mkS_lagrange_localspl_quad		(degree,&rwork[0],&nshapes);
  mkS_canonic_quad(degree,
		   &nshapes,
		   &rwork[2*nshapes],
		   &nshapes,
		   &rwork[0],
		   &nshapes,
		   NULL,
		   NULL,
		   err_);
  mkS_dy_canonic_quad(degree,
		      n,
		      r,
		      roff_,
		      p,
		      poff_,
		      NULL,
		      NULL,
		      err_);
  nslapack_dgesv(&nshapes,
		 n,
		 &rwork[2*nshapes],
		 &nshapes,perm,r,roff_,&info);

}

 void mkS_dx_lagrange_quad(cst_mpsint degree,
				       cst_mpsint n,
				       mpsreal r,
				       cst_mpsint roff_,
				       cst_mpsreal p,
				       cst_mpsint  poff_,
				       mpsreal rwork,
				       cst_mpsint rwork_n,
				       mpsint err_)
{

  nsINT info,nshapes = (degree[0]+1)*(degree[0]+1);
  nsINT perm[128];
  mkS_lagrange_localspl_quad	(degree,&rwork[0],&nshapes);
  mkS_canonic_quad	(degree,&nshapes,&rwork[2*nshapes],&nshapes,&rwork[0],&nshapes,NULL,NULL,err_);
  mkS_dx_canonic_quad	(degree,n,r,roff_,p,poff_,NULL,NULL,err_);
  nslapack_dgesv				(&nshapes,n,&rwork[2*nshapes],&nshapes,perm,r,roff_,&info);
}


void mkS_lagrange_quad(cst_mpsint  degree,
				   cst_mpsint n,
				   mpsreal  r,
				   cst_mpsint roff_,
				   cst_mpsreal  p,
				   cst_mpsint poff_,
				   mpsreal  rwork,
				   cst_mpsint rwork_n,
				   mpsint err_)
{
  nsINT info,nshapes = (degree[0]+1)*(degree[0]+1);
  nsINT perm[128];  
  mkS_lagrange_localspl_quad(degree,&rwork[0],&nshapes);
  mkS_canonic_quad(degree,&nshapes,&rwork[2*nshapes],&nshapes,&rwork[0],&nshapes,NULL,NULL,err_);
  mkS_canonic_quad(degree,n,r,roff_,p,poff_,NULL,NULL,err_);
  nslapack_dgesv(&nshapes,n,&rwork[2*nshapes],&nshapes,perm,r,roff_,&info);  
}
