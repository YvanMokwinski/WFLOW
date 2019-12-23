#include "ns_config_lapack.h"

#if 0
void mkS_lagrange_localspl_interval( cst_pI 	degree_,
				     pR		coo_,
				     cst_pI		off_)
{
  __mmkS_lagrange_localspl_interval(degree_,coo_,off_);
}
#endif

void mkS_dx_lagrange_interval(cst_pI  	degree,
			      cst_pI 	n,
			      pR 	r,
			      cst_pI 	roff_,
			      cst_pR 	p,
			      cst_pI 	poff_,
			      pR 	rwork,
			      cst_pI 	rwork_n,
			      pI 	err_)
{
  err_[0] = (I)0;
  const I k = degree[0];
  switch(k)
    {
    case 1:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)-0.5);
	      r[roff_[0]*i+1] = ((R)0.5);
	    } }
	break;
      }
    case 2:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = p[poff_[0]*i] - ((R)0.5);
	      r[roff_[0]*i+1] = p[poff_[0]*i] + ((R)0.5);
	      r[roff_[0]*i+2] = ((R)-2.0)*p[poff_[0]*i];
	    } }
	break;
      }
    case 0 :
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)0.0);
	    } }
	break;
      }
    case 3:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ( ( ((R)-27.0)*p[poff_[0]*i] + ((R)4.5) )*p[poff_[0]*i] + ((R)1.0) )*((R)0.0625);
	      r[roff_[0]*i+1] = ( ( ((R)27.0)*p[poff_[0]*i] + ((R)4.5) )*p[poff_[0]*i] - ((R)1.0) )*((R)0.0625);
	      r[roff_[0]*i+2] = ( ( ((R)81.0)*p[poff_[0]*i] - ((R)4.5) )*p[poff_[0]*i] - ((R)27.0) )*((R)0.0625);
	      r[roff_[0]*i+3] = ( ( ((R)-81.0)*p[poff_[0]*i] - ((R)4.5) )*p[poff_[0]*i] + ((R)27.0) )*((R)0.0625);
	    } }
	break;
      }
    default:
      {
	const I nshapes = degree[0]+1;
	mkS_lagrange_localspl_interval(degree,
				       &rwork[0],
				       &negal1);
	mkS_canonic_interval	(degree,
				 &nshapes,
				 &rwork[nshapes],
				 &nshapes,
				 &rwork[0],
				 &negal1,
				 NULL,
				 NULL,
				 err_);
	mkS_dx_canonic_interval	(degree,
				 n,
				 r,
				 roff_,
				 p,
				 poff_,
				 NULL,
				 NULL,
				 err_);
	I lapack_info=(I)0,perm[64];
#ifndef NDEBUG
	if (nshapes>64)
	  fprintf(stderr,"hey gros con ! degree = "ifmt" ca te pose pas un ptit probleme ???",degree[0]);
#endif
	dgesv(&nshapes,n,&rwork[nshapes],&nshapes,perm,r,roff_,&lapack_info);  
	break;
      }
    }
}

void mkS_lagrange_interval(cst_pI  	degree,
			   cst_pI 	n,
			   pR 		r,
			   cst_pI 	roff_,
			   cst_pR 	p,
			   cst_pI 	poff_,
			   pR 		rwork,
			   cst_pI 	rwork_n,
			   pI 		err_)
{
  err_[0] 	= (I)0;
  const I k = degree[0];
  switch(k)
    {
    case 1:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ( ((R)1.0) - p[poff_[0]*i] )/((R)2.0);
	      r[roff_[0]*i+1] = ( ((R)1.0) + p[poff_[0]*i] )/((R)2.0);
	    } }

	break;
      }
    case 2:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = p[poff_[0]*i]*(p[poff_[0]*i]-((R)1.0))*((R)0.5);
	      r[roff_[0]*i+1] = p[poff_[0]*i]*(p[poff_[0]*i]+((R)1.0))*((R)0.5);
	      r[roff_[0]*i+2] = ((R)1.0)-p[poff_[0]*i]*p[poff_[0]*i];
	    } }
	break;
      }
    case 0 :
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = ((R)1.0);
	    } }
	break;
      }
    case 3:
      {
	{ I i;
	  for (i=0;i<n[0];++i)
	    {
	      r[roff_[0]*i+0] = (((R)1.0)-p[poff_[0]*i])*(((R)3.0)*p[poff_[0]*i]-((R)1.0))*(((R)3.0)*p[poff_[0]*i]+((R)1.0))*((R)0.0625);
	      r[roff_[0]*i+1] = (((R)1.0)+p[poff_[0]*i])*(((R)3.0)*p[poff_[0]*i]-((R)1.0))*(((R)3.0)*p[poff_[0]*i]+((R)1.0))*((R)0.0625);
	      r[roff_[0]*i+2] = (p[poff_[0]*i]-((R)1.0))*(((R)3.0)*p[poff_[0]*i]-((R)1.0))*(p[poff_[0]*i]+((R)1.0))*((R)0.5625);
	      r[roff_[0]*i+3] = (((R)1.0)-p[poff_[0]*i])*(((R)3.0)*p[poff_[0]*i]+((R)1.0))*(p[poff_[0]*i]+((R)1.0))*((R)0.5625);
	    } }
	break;
      }
    default:
      {
	const I nshapes = degree[0]+1;
	mkS_lagrange_localspl_interval	(degree,
					 &rwork[0],
					 &negal1);
	mkS_canonic_interval		(degree,
					 &nshapes,
					 &rwork[nshapes],
					 &nshapes,
					 &rwork[0],
					 &negal1,
					 NULL,
					 NULL,
					 err_);
	mkS_canonic_interval		(degree,
					 n,
					 r,
					 roff_,
					 p,
					 poff_,
					 NULL,
					 NULL,
					 err_);
	I lapack_info=(I)0,perm[64]; /* celui qui utilise du degre > 63 on le gifle ... */
#ifndef NDEBUG
	if (nshapes>64)
	  fprintf(stderr,"hey gros con ! degree = "ifmt" ca te pose pas un ptit probleme ???\n",degree[0]);
	exit(1);
#endif
	dgesv(&nshapes,n,&rwork[nshapes],&nshapes,perm,r,roff_,&lapack_info);  
	break;
      }
    }
}
