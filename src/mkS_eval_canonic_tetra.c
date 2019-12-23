#include "mkSYS.h"
#include "mkS.h"
#include "mok_m.h"


void mkS_canonic_tetra(cst_mpsint 	degree,
			     cst_mpsint 	n,
			     mpsreal 		r,
			     cst_mpsint 	roff_,
			     cst_mpsreal 	p,
			     cst_mpsint 	poff_,
			     mpsreal 		rwork,
			     cst_mpsint 	rwork_n,
			     mpsint 		err_)
{
  nsINT i;
  err_[0] = (nsINT)0;
  if (degree[0]==0)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)1.0;
    }
  else if (degree[0]==1)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)1.0;
      nsblas_dcopy(n,(nsREAL*)&p[0],&__vmps_blas_negal1,&r[1],roff_);
      nsblas_dcopy(n,(nsREAL*)&p[poff_[0]],&__vmps_blas_negal1,&r[2],roff_);
      nsblas_dcopy(n,(nsREAL*)&p[2*poff_[0]],&__vmps_blas_negal1,&r[3],roff_);      
    }
  else if (degree[0]==2)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)1.0;
      nsblas_dcopy(n,(nsREAL*)&p[0],&__vmps_blas_negal1,&r[1],roff_);
      nsblas_dcopy(n,(nsREAL*)&p[poff_[0]],&__vmps_blas_negal1,&r[2],roff_);
      nsblas_dcopy(n,(nsREAL*)&p[2*poff_[0]],&__vmps_blas_negal1,&r[3],roff_);            
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+4] = p[i]*p[i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+5] = p[i]*p[poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+6] = p[i]*p[2*poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+7] = p[poff_[0]+i]*p[poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+8] = p[poff_[0]+i]*p[2*poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+9] = p[2*poff_[0]+i]*p[2*poff_[0]+i];
    }  
}


void mkS_dx_canonic_tetra(cst_mpsint 	degree,
				cst_mpsint 	n,
				mpsreal 		r,
				cst_mpsint 	roff_,
				cst_mpsreal 	p,
				cst_mpsint 	poff_,
				mpsreal 		rwork,
				cst_mpsint 	rwork_n,
				mpsint 		err_)
{
  nsINT i;
  err_[0] = (nsINT)0;
  if (degree[0]==0)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)0.0;
    }
  else if (degree[0]==1)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+1] = (nsREAL)1.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+2] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+3] = (nsREAL)0.0;
    }
  else if (degree[0]==2)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+1] = (nsREAL)1.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+2] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+3] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+4] = p[i]*((nsREAL)2.0);
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+5] = p[poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+6] = p[2*poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+7] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+8] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+9] = (nsREAL)0.0;
    }  
}

void mkS_dy_canonic_tetra(cst_mpsint 	degree,
				cst_mpsint 	n,
				mpsreal 		r,
				cst_mpsint 	roff_,
				cst_mpsreal 	p,
				cst_mpsint 	poff_,
				mpsreal 		rwork,
				cst_mpsint 	rwork_n,
				mpsint 		err_)
{
  nsINT i;
  err_[0] = (nsINT)0;
  if (degree[0]==0)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)0.0;
    }
  else if (degree[0]==1)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+1] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+2] = (nsREAL)1.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+3] = (nsREAL)0.0;
    }
  else if (degree[0]==2)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+1] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+2] = (nsREAL)1.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+3] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+4] = ((nsREAL)0.0);
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+5] = p[i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+6] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+7] = p[poff_[0]+i]*((nsREAL)2.0);
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+8] = p[poff_[0]*2+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+9] = (nsREAL)0.0;
    }  
}

void mkS_dz_canonic_tetra(cst_mpsint 	degree,
				cst_mpsint 	n,
				mpsreal 		r,
				cst_mpsint 	roff_,
				cst_mpsreal 	p,
				cst_mpsint 	poff_,
				mpsreal 		rwork,
				cst_mpsint 	rwork_n,
				mpsint 		err_)
{
  nsINT i;
  err_[0] = (nsINT)0;
  if (degree[0]==0)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)0.0;
    }
  else if (degree[0]==1)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+1] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+2] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+3] = (nsREAL)1.0;
    }
  else if (degree[0]==2)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+1] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+2] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+3] = (nsREAL)1.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+4] = ((nsREAL)0.0);
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+5] = ((nsREAL)0.0);
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+6] = p[i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+7] = ((nsREAL)0.0) ;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+8] = p[poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+9] = p[poff_[0]*2+i]*((nsREAL)2.0);
    }  
}

