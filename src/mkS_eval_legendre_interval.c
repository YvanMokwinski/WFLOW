#include "mkSYS.h"
#include "mkS.h"
#include "mok_m.h"


/*
  CALCUL PAR UNE FORMULE DE RECURRENCE DE BONNET
*/
void mkS_legendre_interval(cst_mkI  	degree,
			   cst_mkI 	n,
			   mkR 		r,
			   cst_mkI 	roff_,
			   cst_mkR 	p,
			   cst_mkI 	poff_,
			   mkR 		rwork,
			   cst_mkI 	rwork_n,
			   mkI 		err_)
{
  err_[0] = (nsINT)NO;
  { const nsINT k = degree[0];
    { nsINT i;
      for (i=0;i<n[0];++i)
	{
	  r[roff_[0]*i+0] = ((nsREAL)1.0);
	  r[roff_[0]*i+1] = p[poff_[0]*i];
	  { nsINT j;
	    for (j=1;j<k;++j)
	      { 
		r[roff_[0]*i+(j+1)] = (p[poff_[0]*i]*((nsREAL)(2*j+1))*r[roff_[0]*i+j] - r[roff_[0]*i+(j-1)]*((nsREAL)j)  )/((nsREAL)(j+1));
	      } }
	} } }
}

void mkS_dx_legendre_interval(cst_mkI  	degree,
			      cst_mkI 	n,
			      mkR 	r,
			      cst_mkI 	roff_,
			      cst_mkR 	p,
			      cst_mkI 	poff_,
			      mkR 	rwork,
			      cst_mkI 	rwork_n,
			      mkI 	err_)
{
  err_[0] = (nsINT)NO;
  { const nsINT k = degree[0];
    { nsINT i;
      for (i=0;i<n[0];++i)
	{
	  r[roff_[0]*i+0] = ((nsREAL)0.0);
	  r[roff_[0]*i+1] = ((nsREAL)1.0);
	  nsREAL pp0 	  = ((nsREAL)1.0);
	  nsREAL pp1 	  = p[poff_[0]*i];
	  nsREAL pp2 	  = p[poff_[0]*i];
	  { nsINT j;
	    for (j=1;j<k;++j)
	      { 
		r[roff_[0]*i+(j+1)] 	= (p[poff_[0]*i]*((nsREAL)(2*j+1))*r[roff_[0]*i+j] - r[roff_[0]*i+(j-1)]*((nsREAL)j)  )/((nsREAL)(j+1)) + pp2 * (((nsREAL)(2*j+1))) / (((nsREAL)(j+1)));
		{ const nsREAL tmp 	= pp2;
		  pp2 			= (p[poff_[0]*i]*((nsREAL)(2*j+1))*pp1 - pp0*((nsREAL)j)  )/((nsREAL)(j+1));
		  pp0 			= pp1;
		  pp1 			= tmp; }
	      } }
	} } }
}
