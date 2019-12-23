#include "mkSYS.h"
#include "mkS.h"
#include "mok_m.h"

void mkS_canonic_quad(      cst_mpsint 	degree,
			    cst_mpsint 	n,
			    mpsreal 	r,
			    cst_mpsint 	roff_,
			    cst_mpsreal p,
			    cst_mpsint 	poff_,
			    mpsreal 	rwork,
			    cst_mpsint 	rwork_n,
			    mpsint 	err_)
{
  nsINT negal1 = (nsINT)1;
  nsINT i,j,k,s;
  err_[0] = 0;
  for (i=0;i<n[0];++i)
    r[i*roff_[0]] = (nsREAL)1.0;
  if (degree[0]>0)
    {
      nsblas_dcopy(n,(nsREAL*)&p[0],&negal1,&r[1],roff_);
      nsblas_dcopy(n,(nsREAL*)&p[poff_[0]],&negal1,&r[2],roff_);
      for (i=0;i<n[0];++i)
	r[3+i*roff_[0]] = r[1+i*roff_[0]]*r[2+i*roff_[0]];   
      s = 4;
      for (k=2;k<=degree[0];++k)
	{
	  for (j=0;j<=k;++j,++s)
	    for (i=0;i<n[0];++i)
	      r[s+i*roff_[0]] = nsPOW(p[i],(nsREAL)k) * nsPOW(p[poff_[0]+i],(nsREAL)j);	  
	  for (j=1;j<=k;++j,++s)
	    for (i=0;i<n[0];++i)
	      r[s+i*roff_[0]] = nsPOW(p[i],(nsREAL)(k-j)) * nsPOW(p[poff_[0]+i],(nsREAL)k);	 	  
	}
#if __mk_debug__
      if (s!=(degree[0]+1)*(degree[0]+1))
	mk_msg("mkS_canonic_quad failed");
#endif
    }
}

void mkS_dx_canonic_quad(cst_mpsint 	degree,
			 cst_mpsint 	n,
			 mpsreal 	r,
			 cst_mpsint 	roff_,
			 cst_mpsreal 	p,
			 cst_mpsint 	poff_,
			 mpsreal 	rwork,
			 cst_mpsint 	rwork_n,
			 mpsint 	err_)
{
  nsINT i,j,k,s;
  err_[0] = 0;
  for (i=0;i<n[0];++i)
    r[i*roff_[0]] = (nsREAL)0.0;
  if (degree[0]>0)
    {
      for (i=0;i<n[0];++i)
	r[1+i*roff_[0]] = (nsREAL)1.0;
      for (i=0;i<n[0];++i)
	r[2+i*roff_[0]] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[3+i*roff_[0]] = p[poff_[0]+i];   
      s = 4;
      for (k=2;k<=degree[0];++k)
	{
	  for (j=0;j<=k;++j,++s)
	    for (i=0;i<n[0];++i)
	      r[s+i*roff_[0]] = ((nsREAL)k) * nsPOW(p[i],(nsREAL)(k-1)) * nsPOW(p[poff_[0]+i],(nsREAL)j);
	  
	  for (j=1;j<=k;++j,++s)
	    for (i=0;i<n[0];++i)
	      r[s+i*roff_[0]] = ((nsREAL)(k-j))*nsPOW(p[i],(nsREAL)(k-j-1)) * nsPOW(p[poff_[0]+i],(nsREAL)k);	 	  
	}
#if __mk_debug__
      if (s!=(degree[0]+1)*(degree[0]+1))
	mk_msg("mkS_dx_canonic_quad failed");
#endif
    }
  
}

void mkS_dy_canonic_quad(cst_mpsint 	degree,
			 cst_mpsint 	n,
			 mpsreal 	r,
			 cst_mpsint 	roff_,
			 cst_mpsreal 	p,
			 cst_mpsint 	poff_,
			 mpsreal 	rwork,
			 cst_mpsint 	rwork_n,
			 mpsint 	err_)
{
  nsINT i,j,k,s;
  err_[0] = 0;
  for (i=0;i<n[0];++i)
    r[i*roff_[0]] = (nsREAL)0.0;
  if (degree[0]>0)
    {
      for (i=0;i<n[0];++i)
	r[1+i*roff_[0]] = (nsREAL)0.0;
      for (i=0;i<n[0];++i)
	r[2+i*roff_[0]] = (nsREAL)1.0;
      for (i=0;i<n[0];++i)
	r[3+i*roff_[0]] = p[i];   
      s = 4;
      for (k=2;k<=degree[0];++k)
	{
	  for (j=0;j<=k;++j,++s)
	    for (i=0;i<n[0];++i)
	      r[s+i*roff_[0]] = ((nsREAL)j) * nsPOW(p[i],(nsREAL)k) * nsPOW(p[poff_[0]+i],(nsREAL)(j-1));
	  
	  for (j=1;j<=k;++j,++s)
	    for (i=0;i<n[0];++i)
	      r[s+i*roff_[0]] = ((nsREAL)k) * nsPOW(p[i],(nsREAL)(k-j)) * nsPOW(p[poff_[0]+i],(nsREAL)(k-1));	 	  
	}
#if __mk_debug__
      if (s!=(degree[0]+1)*(degree[0]+1))
	mk_msg("mkS_dy_canonic_quad failed");
#endif
    }
}




