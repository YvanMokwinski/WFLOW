#include "mkS.h"
#include <math.h>
#if 0
#include "mkSYS.h"

#include "mok_m.h"
#endif
void mkS_canonic_tria(cst_pI 	degree,
			    cst_pI 	n,
			    pR 		r,
			    cst_pI 	roff_,
			    cst_pR 	p,
			    cst_pI 	poff_,
			    pR 		rwork,
			    cst_pI 	rwork_n,
			    pI 		err_)
{
  I i,q,k,j;
  err_[0] = (I)0;
  if (degree[0]==0)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (R)1.0;
    }
  else if (degree[0]==1)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];
	  r[i*roff_[0]] = (R)1.0;
	  r[i*roff_[0]+1] = x;
	  r[i*roff_[0]+2] = y;
	}
    }
  else if (degree[0]==2)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];
	  r[i*roff_[0]] = (R)1.0;
	  r[i*roff_[0]+1] = x;
	  r[i*roff_[0]+2] = y;
	  r[i*roff_[0]+3] = x*x;
	  r[i*roff_[0]+4] = y*x;
	  r[i*roff_[0]+5] = y*y;
	}
    }
  else if (degree[0]==3)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];
	  r[i*roff_[0]] = (R)1.0;
	  r[i*roff_[0]+1] = x;
	  r[i*roff_[0]+2] = y;
	  r[i*roff_[0]+3] = x*x;
	  r[i*roff_[0]+4] = y*x;
	  r[i*roff_[0]+5] = y*y;
	  r[i*roff_[0]+6] = r[i*roff_[0]+3]*x;
	  r[i*roff_[0]+7] = r[i*roff_[0]+3]*y;
	  r[i*roff_[0]+8] = r[i*roff_[0]+5]*x;
	  r[i*roff_[0]+9] = r[i*roff_[0]+5]*y;
	}
    }
  else if (degree[0]==4)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x 	= p[i];
	  const R y 	= p[poff_[0]+i];
	  r[i*roff_[0]] 	= (R)1.0;
	  r[i*roff_[0]+1] 	= x;
	  r[i*roff_[0]+2] 	= y;
	  r[i*roff_[0]+3] 	= x*x;
	  r[i*roff_[0]+4] 	= y*x;
	  r[i*roff_[0]+5] 	= y*y;
	  r[i*roff_[0]+6] 	= r[i*roff_[0]+3]*x;
	  r[i*roff_[0]+7] 	= r[i*roff_[0]+3]*y;
	  r[i*roff_[0]+8] 	= r[i*roff_[0]+5]*x;
	  r[i*roff_[0]+9] 	= r[i*roff_[0]+5]*y;
	  r[i*roff_[0]+10] 	= r[i*roff_[0]+6]*x;
	  r[i*roff_[0]+11] 	= r[i*roff_[0]+6]*y;
	  r[i*roff_[0]+12] 	= r[i*roff_[0]+3]*r[i*roff_[0]+5];
	  r[i*roff_[0]+13] 	= r[i*roff_[0]+9]*x;
	  r[i*roff_[0]+14] 	= r[i*roff_[0]+9]*y;
	}
    }
  else
    {
      for (i=0;i<n[0];++i)
	{
	  const R x 	= p[i];
	  const R y 	= p[poff_[0]+i];
	  r[i*roff_[0]] = (R)1.0;
	  r[i*roff_[0]+1] 	= x;
	  r[i*roff_[0]+2] 	= y;
	  r[i*roff_[0]+3] 	= x*x;
	  r[i*roff_[0]+4] 	= y*x;
	  r[i*roff_[0]+5] 	= y*y;
	  r[i*roff_[0]+6] 	= r[i*roff_[0]+3]*x;
	  r[i*roff_[0]+7] 	= r[i*roff_[0]+3]*y;
	  r[i*roff_[0]+8] 	= r[i*roff_[0]+5]*x;
	  r[i*roff_[0]+9] 	= r[i*roff_[0]+5]*y;
	  r[i*roff_[0]+10] 	= r[i*roff_[0]+6]*x;
	  r[i*roff_[0]+11] 	= r[i*roff_[0]+6]*y;
	  r[i*roff_[0]+12] 	= r[i*roff_[0]+3]*r[i*roff_[0]+5];
	  r[i*roff_[0]+13] 	= r[i*roff_[0]+9]*x;
	  r[i*roff_[0]+14] 	= r[i*roff_[0]+9]*y;
	}
      for (k=5;k<=degree[0];++k)
	for (q = (k*(k+1))/((I)2),j=0;j<=k;++j)
	  for (i=0;i<n[0];++i)
	    r[i*roff_[0]+q+j] = nsPOW(p[i],((R)k-j)) * nsPOW(p[poff_[0]+i],((R)j));    
    }
}


void mkS_dy_canonic_tria(cst_pI 	degree,
			    cst_pI 	n,
			    pR 		r,
			    cst_pI 	roff_,
			    cst_pR 	p,
			    cst_pI 	poff_,
			    pR 		rwork,
			    cst_pI 	rwork_n,
			    pI 		err_)
{
  I i,q,k,j;
  err_[0] = (I)0;
  if (degree[0]==0)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (R)0.0;
    }
  else if (degree[0]==1)
    {
      for (i=0;i<n[0];++i)
	{
	  r[i*roff_[0]] = (R)0.0;
	  r[i*roff_[0]+1] = (R)0.0;
	  r[i*roff_[0]+2] = (R)1.0;
	}
    }
  else if (degree[0]==2)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];
	  r[i*roff_[0]] = (R)0.0;
	  r[i*roff_[0]+1] = (R)0.0;
	  r[i*roff_[0]+2] = (R)1.0;
	  r[i*roff_[0]+3] = (R)0.0;
	  r[i*roff_[0]+4] = x;
	  r[i*roff_[0]+5] = y*((R)2.0);
	}
    }
  else if (degree[0]==3)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];
	  r[i*roff_[0]] = (R)0.0;
	  r[i*roff_[0]+1] = (R)0.0;
	  r[i*roff_[0]+2] = (R)1.0;
	  r[i*roff_[0]+3] = (R)0.0;
	  r[i*roff_[0]+4] = x;
	  r[i*roff_[0]+5] = y*((R)2.0);
	  r[i*roff_[0]+6] = ((R)0.0);
	  r[i*roff_[0]+7] = x*x;
	  r[i*roff_[0]+8] = ((R)2.0)*x*y;
	  r[i*roff_[0]+9] = ((R)3.0)*y*y;
	}
    }
  else if (degree[0]==4)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];

	  r[i*roff_[0]] = (R)0.0;
	  r[i*roff_[0]+1] = (R)0.0;
	  r[i*roff_[0]+2] = (R)1.0;
	  r[i*roff_[0]+3] = (R)0.0;
	  r[i*roff_[0]+4] = x;
	  r[i*roff_[0]+5] = y*((R)2.0);
	  r[i*roff_[0]+6] = ((R)0.0);
	  r[i*roff_[0]+7] = x*x;
	  r[i*roff_[0]+8] = ((R)2.0)*x*y;
	  r[i*roff_[0]+9] = ((R)3.0)*y*y;

	  r[i*roff_[0]+10] = ((R)0.0);
	  r[i*roff_[0]+11] = ((R)1.0)*x*x*x;
	  r[i*roff_[0]+12] = ((R)2.0)*y*x*x;
	  r[i*roff_[0]+13] = ((R)3.0)*y*y*x;
	  r[i*roff_[0]+14] = ((R)4.0)*y*y*y;
	}
    }
  else
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];

	  r[i*roff_[0]] = (R)0.0;
	  r[i*roff_[0]+1] = (R)0.0;
	  r[i*roff_[0]+2] = (R)1.0;
	  r[i*roff_[0]+3] = (R)0.0;
	  r[i*roff_[0]+4] = x;
	  r[i*roff_[0]+5] = y*((R)2.0);
	  r[i*roff_[0]+6] = ((R)0.0);
	  r[i*roff_[0]+7] = x*x;
	  r[i*roff_[0]+8] = ((R)2.0)*x*y;
	  r[i*roff_[0]+9] = ((R)3.0)*y*y;

	  r[i*roff_[0]+10] = ((R)0.0);
	  r[i*roff_[0]+11] = ((R)1.0)*x*x*x;
	  r[i*roff_[0]+12] = ((R)2.0)*y*x*x;
	  r[i*roff_[0]+13] = ((R)3.0)*y*y*x;
	  r[i*roff_[0]+14] = ((R)4.0)*y*y*y;
	}
      for (k=5;k<=degree[0];++k)
	{ 
	  q = (k*(k+1))/((I)2);
	  for (i=0;i<n[0];++i)
	    r[i*roff_[0]+q+0] = (R)0.0;
	  for (j=1;j<=k;++j)
	    for (i=0;i<n[0];++i)
	      r[i*roff_[0]+q+j] = ( ((R)j) * nsPOW(p[poff_[0]+i],((R)j-1)) ) * nsPOW(p[i],((R)(k-j)));
	}
    }
}


void mkS_dx_canonic_tria(cst_pI 	degree,
			    cst_pI 	n,
			    pR 		r,
			    cst_pI 	roff_,
			    cst_pR 	p,
			    cst_pI 	poff_,
			    pR 		rwork,
			    cst_pI 	rwork_n,
			    pI 		err_)
{
  I i,q,k,j;
  err_[0] = (I)0;
  if (degree[0]==0)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (R)0.0;
    }
  else if (degree[0]==1)
    {
      for (i=0;i<n[0];++i)
	{
	  r[i*roff_[0]] 	= (R)0.0;
	  r[i*roff_[0]+1] 	= (R)1.0;
	  r[i*roff_[0]+2] 	= (R)0.0;
	}
    }
  else if (degree[0]==2)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];
	  r[i*roff_[0]] = (R)0.0;
	  r[i*roff_[0]+1] = (R)1.0;
	  r[i*roff_[0]+2] = (R)0.0;
	  r[i*roff_[0]+3] = ((R)2.0)*x;
	  r[i*roff_[0]+4] = y;
	  r[i*roff_[0]+5] = (R)0.0;
	}
    }
  else if (degree[0]==3)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];
	  r[i*roff_[0]] = (R)0.0;
	  r[i*roff_[0]+1] = (R)1.0;
	  r[i*roff_[0]+2] = (R)0.0;
	  r[i*roff_[0]+3] = ((R)2.0)*x;
	  r[i*roff_[0]+4] = y;
	  r[i*roff_[0]+5] = (R)0.0;

	  r[i*roff_[0]+6] = ((R)3.0)*x*x;
	  r[i*roff_[0]+7] = ((R)2.0)*x*y;
	  r[i*roff_[0]+8] = y*y;
	  r[i*roff_[0]+9] = ((R)0.0);
	}
    }
  else if (degree[0]==4)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];
	  r[i*roff_[0]] = (R)0.0;
	  r[i*roff_[0]+1] = (R)1.0;
	  r[i*roff_[0]+2] = (R)0.0;
	  r[i*roff_[0]+3] = ((R)2.0)*x;
	  r[i*roff_[0]+4] = y;
	  r[i*roff_[0]+5] = (R)0.0;

	  r[i*roff_[0]+6] = ((R)3.0)*x*x;
	  r[i*roff_[0]+7] = ((R)2.0)*x*y;
	  r[i*roff_[0]+8] = y*y;
	  r[i*roff_[0]+9] = ((R)0.0);

	  r[i*roff_[0]+10] = ((R)4.0)*x*x*x;
	  r[i*roff_[0]+11] = ((R)3.0)*x*x*y;
	  r[i*roff_[0]+12] = ((R)2.0)*y*y*x;
	  r[i*roff_[0]+13] = ((R)1.0)*y*y*y;
	  r[i*roff_[0]+14] = ((R)0.0);
	}
    }
  else
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];
	  r[i*roff_[0]] = (R)0.0;
	  r[i*roff_[0]+1] = (R)1.0;
	  r[i*roff_[0]+2] = (R)0.0;

	  r[i*roff_[0]+3] = ((R)2.0)*x;
	  r[i*roff_[0]+4] = y;
	  r[i*roff_[0]+5] = (R)0.0;

	  r[i*roff_[0]+6] = ((R)3.0)*x*x;
	  r[i*roff_[0]+7] = ((R)2.0)*x*y;
	  r[i*roff_[0]+8] = y*y;
	  r[i*roff_[0]+9] = ((R)0.0);

	  r[i*roff_[0]+10] = ((R)4.0)*x*x*x;
	  r[i*roff_[0]+11] = ((R)3.0)*x*x*y;
	  r[i*roff_[0]+12] = ((R)2.0)*y*y*x;
	  r[i*roff_[0]+13] = ((R)1.0)*y*y*y;
	  r[i*roff_[0]+14] = ((R)0.0);
	}
      for (k=5;k<=degree[0];++k)
	{  
	  q = (k*(k+1))/((I)2);
	  for (j=0;j<k;++j)
	    for (i=0;i<n[0];++i)
	      r[i*roff_[0]+q+j] = ( ((R)(k-j)) * nsPOW(p[i],((R)(k-j-1))) ) * nsPOW(p[poff_[0]+i],((R)j));
	  for (i=0;i<n[0];++i)
	    r[i*roff_[0]+q+j] = ((R)0.0);
	}
    }
}


#if 0
void mkS_dx_canonic_tria(cst_pI 	degree,
			       cst_pI 	n,
			       pR		r,
			       cst_pI 	roff_,
			       cst_pR 	p,
			       cst_pI 	poff_,
			       pR		rwork,
			       cst_pI 	rwork_n,
			       pI 		err_)
{  I i,j,k,q ;
  err_ [0] = (I)0;
  for (i=0;i<n[0];++i)
    r[i*roff_[0]] = (R)0.0;
  if (degree[0]>0)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+1] = (R)1.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+2] = (R)0.0;
      for (k=2;k<=degree[0];++k)
	{  
	  q = (k*(k+1))/((I)2);
	  for (j=0;j<k;++j)
	    for (i=0;i<n[0];++i)
	      r[i*roff_[0]+q+j] = ( ((R)(k-j)) * nsPOW(p[i],((R)(k-j-1))) ) * nsPOW(p[poff_[0]+i],((R)j));
	  for (i=0;i<n[0];++i)
	    r[i*roff_[0]+q+j] = 0.0;
	}
    }
}

void mkS_dy_canonic_tria(cst_pI 	degree,
			       cst_pI 	n,
			       pR	r,
			       cst_pI 	roff_,
			       cst_pR 	p,
			       cst_pI 	poff_,
			       R*	rwork,
			       cst_pI 	rwork_n,
			       pI 	err_)
{
  err_[0] = (I)0;
  I i,j,k,q;
  for (i=0;i<n[0];++i)
    r[i*roff_[0]] = (R)0.0;
  if (degree[0]>0)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+1] = (R)0.0;
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+2] = (R)1.0;
      for (k=2;k<=degree[0];++k)
	{  
	  q = (k*(k+1))/((I)2);
	  for (j=0,i=0;i<n[0];++i)
	    r[i*roff_[0]+q+j] = (R)0.0;
	  for (j=1;j<=k;++j)
	    for (i=0;i<n[0];++i)
	      r[i*roff_[0]+q+j] = ( ((R)j) * nsPOW(p[poff_[0]+i],((R)(j-1))) ) * nsPOW(p[i],((R)(k-j)));
	}
    }
}


#endif
