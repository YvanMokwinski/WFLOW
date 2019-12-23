
#include "mkS.h"
#include "Blas.h"
#if 0
#include "mkSYS.h"

#include "mok_m.h"
#include "mkS_eval.h"
void mkS_solve_meilleur()
{
  R cooelm[6],P[_N*dimension],R[_N],ih_c[_N*_N],cc_tmp[_N*_N],gl_c[_N*2],C[_N],laprwork[_N*4],matp[_N*_N],rhsp[_N*dimension];
  static const char FACT[1] = {'E'};
  static char EQUED[1]={'B'};
  R RCOND[1],FERR[1],BERR[1],AF[_N*_N],MINRCOND=1.0e+30,MAXRCOND=0.0,MINFERR=1.0e+30,MAXFERR=0.0,MINBERR=1.0e+30,MAXBERR=0.0;
  dgesv(&nshapes,n,&rwork[2*nshapes],&nshapes,perm,r,roff_,&info);
  dgesvx_ (FACT,
	   transN,
	   &nshapes,
	   (NOT zz_[0])?(&n1):(&dim),
	   matp,
	   &nshapes,
	   rwork[2*nshapes],
	   &nshapes,
	   perm,
	   EQUED,
	   R,
	   C,
	   rhsp,
	   &N,
	   P,
	   &N,
	   RCOND,
	   FERR,
	   BERR,
	   laprwork,
	   lapiwork,
	   info_);  
  if (RCOND[0]<MINRCOND)
    {
      MINRCOND=RCOND[0];
    }
  if (RCOND[0]>MAXRCOND)
    {
      MAXRCOND=RCOND[0];
    }  
  if (FERR[0]<MINFERR)
    {
      MINFERR=FERR[0];
    }
  if (FERR[0]>MAXFERR)
    {
      MAXFERR=FERR[0];
    }
  if (BERR[0]<MINBERR)
    {
      MINBERR=BERR[0];
    }
  if (BERR[0]>MAXBERR)
    {
      MAXBERR=BERR[0];
    }
}
#endif


#define resoluteur_system dgesv

#if 0

/* CA FONCTIONNE */
#include "mkFACTORIZATION.h"

void ssolve(cst_pI n_,
	    cst_pI nrhs,
	    pR A_,
	    cst_pI Aoff_,
	    pI perm,
	    pR rhs_,
	    cst_pI rhsoff,
	    cst_pI info_)
{
  emk_err err_[1]={__emk_err_no};
  static R rwork[32000];
  static I 	ipar[64];
  static R fpar[64],sol[64]; 
  { I i;
    for (i=0;i<nrhs[0];++i)
      {
	I its=0;
	memset(ipar,0,64*sizeof(I));
	memset(fpar,0,64*sizeof(R));
	ipar[2-1] 		= (I)0; 
	ipar[3-1] 		= (I)1; /* stopping criteria */
	ipar[5-1] 		= (I)20;
	ipar[6-1] 		= (I)1000;
	fpar[0] 		= (R)1.0e-10;
	fpar[1] 		= (R)2.22e-16;
	ipar[1-1] 		= (I)0;  
	nsFLAG we_continue = 1;
	while (we_continue)
	  {
	    { I n = n_[0];
	      mok_use_sparsekit_gmres(&n,
				      &rhs_[rhsoff[0]*i],
				      sol,
				      ipar,
				      fpar,
				      rwork); }
	    if (ipar[7-1] != its)
	      {
		its = ipar[7-1];
	      }
	    int kk = ipar[0];
	    switch (kk)
	      {
	      case 1:
		{
		  Blas_dgemv(__vmps_blas_transN,n_,n_,&__vmps_blas_regal1,A_,Aoff_,&rwork[ipar[8-1]-1],&__vmps_blas_negal1,&__vmps_blas_regal0,&rwork[ipar[9-1]-1],&__vmps_blas_negal1);
		  break;
		}
	      case 2:
		{
		  break;
		}
	      case 3:
	      case 5:	  
		{	  
		  break;
		}
	      case 4:
	      case 6:	  
		{
		  break;
		}
	      case 0 :
		{
		  we_continue = (nsFLAG)NO;
		  break;
		}
	      case -1 :
	      case -2:
	      case -3:
		{
		  {
		    err_[0] 	= __emk_err_extern_lib;
		    we_continue = (nsFLAG)NO;
		    we_continue = (nsFLAG)NO;
		    break;
		  }
		}
	      default:
		{
		  err_[0] 	= __emk_err_switch;
		  we_continue = (nsFLAG)NO;
		  break;
		}	  
	      }  
	    if (err_[0])
	      break;
	  }
	Blas_dcopy(n_,sol,&__vmps_blas_negal1,&rhs_[rhsoff[0]*i],&__vmps_blas_negal1);
      } }	
}
#endif




#if __mk_debug__
void mkS_lagrange_localspl_tria( cst_pI 	degree,
				 pR		coo_,
				 cst_pI	off_)
{
  __mmkS_lagrange_localspl_tria(degree,coo_,off_);
}
#endif

void mkS_dy_lagrange_tria(cst_pI degree,
			  cst_pI n,
			  pR r,
			  cst_pI roff_,
			  cst_pR p,
			  cst_pI poff_,
			  pR rwork,
			  cst_pI rwork_n,
			  pI err_)
{
  err_[0] = (I)0;
  if (degree[0]==1)
    {
      { I i;
	for (i=0;i<n[0];++i)
	  {
	    r[roff_[0]*i+1] = ((R)0.0);
	    r[roff_[0]*i+2] = ((R)1.0);
	    r[roff_[0]*i+0] = ((R)-1.0);
	  } }
    }
  else if (degree[0]==0)
    {
      { I i;
	for (i=0;i<n[0];++i)
	  {
	    r[roff_[0]*i+0] = ((R)0.0);
	  } }
    }
  else if (degree[0]==2)
    {
      { I i;
	for (i=0;i<n[0];++i)
	  {
	    const R x   	= p[i];
	    const R y   	= p[poff_[0]+i];
	    const R lam 	= ((R)1.0)-(x+y);
	    r[roff_[0]*i+0] = ((R)1.0)-((R)4.0)*lam;
	    r[roff_[0]*i+1] = (R)0.0;
	    r[roff_[0]*i+2] = ((R)-1.0)+((R)4.0)*y;
	    r[roff_[0]*i+3] = ((R)-4.0)*x;
	    r[roff_[0]*i+4] = ((R)4.0)*x;
	    r[roff_[0]*i+5] = ((R)4.0)*(lam-y);
	  } }
    }
  else if (degree[0]==3)
    {
      { I i;
	for (i=0;i<n[0];++i)
	  {
	    const R x   = p[i];
	    const R y   = p[poff_[0]+i];
	    const R lam = ((R)1.0)-(x+y);

	    r[roff_[0]*i+0]= ((R)0.5)*(((-((R)2.0))+((R)18.0)*lam)-((R)27.0)*lam*lam);
	    r[roff_[0]*i+1]= (R)0.0;
	    r[roff_[0]*i+2]= ((R)0.5)*((((R)2.0)-((R)18.0)*y)+((R)27.0)*y*y);
	    r[roff_[0]*i+3]= ((R)0.5)*(((R)-9.0))*x*((-((R)1.0))+((R)6.0)*lam);
	    r[roff_[0]*i+4]= ((R)0.5)*(((R)-9.0))*x*((-((R)1.0))+((R)3.0)*x);
	    r[roff_[0]*i+5]= ((R)0.5)*((R)9.0)*x*((-((R)1.0))+((R)3.0)*x);
	    r[roff_[0]*i+6]= ((R)0.5)*((R)9.0)*x*((-((R)1.0))+((R)6.0)*y);
	    r[roff_[0]*i+7]= ((R)0.5)*(((R)9.0)*y*((((R)1.0)+((R)6.0)*lam)-((R)3.0)*y)-((R)9.0)*lam);
	    r[roff_[0]*i+8]= ((R)0.5)*(((R)9.0)*lam*(((-((R)1.0))+((R)3.0)*lam)-((R)6.0)*y)+((R)9.0)*y);
	    r[roff_[0]*i+9]= ((R)0.5)*((R)54.0)*x*(lam-y);

	  } }
    }
  else
    {      
      err_[0] = (I)0;
      I info,nshapes = ((degree[0]+1)*(degree[0]+2))/2;
      I perm[nshapes];
      mkS_lagrange_localspl_tria	(degree,&rwork[0],&nshapes);
      mkS_canonic_tria	(degree,
			 &nshapes,
			 &rwork[2*nshapes],
			 &nshapes,
			 &rwork[0],
			 &nshapes,
			 NULL,
			 NULL,
			 err_);  
      mkS_dy_canonic_tria(degree,
			  n,
			  r,
			  roff_,
			  p,
			  poff_,
			  NULL,
			  NULL,
			  err_);  
      dgesv(&nshapes,n,&rwork[2*nshapes],&nshapes,perm,r,roff_,&info);
    }
}


void mkS_dx_lagrange_tria(cst_pI degree,
			  cst_pI n,
			  pR r,
			  cst_pI roff_,
			  cst_pR p,
			  cst_pI poff_,
			  pR rwork,
			  cst_pI rwork_n,
			  pI err_)
{
  err_[0] = (I)0;
  if (degree[0]==1)
    {
      { I i;
	for (i=0;i<n[0];++i)
	  {
	    r[roff_[0]*i+1] = ((R)1.0);
	    r[roff_[0]*i+2] = ((R)0.0);
	    r[roff_[0]*i+0] = ((R)-1.0);
	  } }
    }
  else if (degree[0]==0)
    {
      { I i;
	for (i=0;i<n[0];++i)
	  {
	    r[roff_[0]*i+0] = ((R)0.0);
	  } }
    }
  else if (degree[0]==2)
    {

      { I i;
	for (i=0;i<n[0];++i)
	  {
	    const R x   	= p[i];
	    const R y   	= p[poff_[0]+i];
	    const R lam 	= ((R)1.0)-(x+y);
	    r[roff_[0]*i+0] = ((R)1.0)-((R)4.0)*lam;
	    r[roff_[0]*i+1] = ((R)-1.0)+((R)4.0)*x;
	    r[roff_[0]*i+2] = (R)0.0;
	    r[roff_[0]*i+3] = ((R)4.0)*(lam-x);
	    r[roff_[0]*i+4] = ((R)4.0)*y;  
	    r[roff_[0]*i+5] = ((R)-4.0)*y;
	  } }



    }
  else if (degree[0]==3)
    {           
      { I i;
	for (i=0;i<n[0];++i)
	  {
	    const R x   = p[i];
	    const R y   = p[poff_[0]+i];
	    const R lam = ((R)1.0)-(x+y);
	    
	    r[roff_[0]*i+0]= ((R)0.5)*(((-((R)2.0))+((R)18.0)*lam)-((R)27.0)*lam*lam);
	    r[roff_[0]*i+1]= ((R)0.5)*((((R)2.0)-((R)18.0)*x)+((R)27.0)*x*x);
	    r[roff_[0]*i+2]= (R)0.0;
	    r[roff_[0]*i+3]= ((R)0.5)*(((R)9.0)*lam*(((-((R)1.0))+((R)3.0)*lam)-((R)6.0)*x)+((R)9.0)*x);
	    r[roff_[0]*i+4]= ((R)0.5)*(((R)9.0)*x*((((R)1.0)+((R)6.0)*lam)-((R)3.0)*x)-((R)9.0)*lam);
	    r[roff_[0]*i+5]= ((R)0.5)*((R)9.0)*y*((-((R)1.0))+((R)6.0)*x);
	    r[roff_[0]*i+6]= ((R)0.5)*((R)9.0)*y*((-((R)1.0))+((R)3.0)*y);
	    r[roff_[0]*i+7]= ((R)0.5)*(((R)-9.0))*y*((-((R)1.0))+((R)3.0)*y);
	    r[roff_[0]*i+8]= ((R)0.5)*(((R)-9.0))*y*((-((R)1.0))+((R)6.0)*lam);  
	    r[roff_[0]*i+9]= ((R)0.5)*((R)54.0)*y*(lam-x);
	  } }
    }
  else
    {
      
      err_[0] = (I)0;
      I info,nshapes = ((degree[0]+1)*(degree[0]+2))/2;
      I perm[nshapes];
      mkS_lagrange_localspl_tria	(degree,&rwork[0],&nshapes);
      mkS_canonic_tria	(degree,
			 &nshapes,
			 &rwork[2*nshapes],
			 &nshapes,
			 &rwork[0],
			 &nshapes,
			 NULL,
			 NULL,
			 err_);  
      mkS_dx_canonic_tria(degree,
			  n,
			  r,
			  roff_,
			  p,
			  poff_,
			  NULL,
			  NULL,
			  err_);  
      dgesv(&nshapes,n,&rwork[2*nshapes],&nshapes,perm,r,roff_,&info);
    }
}





void mkS_lagrange_tria(cst_pI  degree,
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
  if (degree[0]==1)
    {
      { I i;
	for (i=0;i<n[0];++i)
	  {
	    r[roff_[0]*i+1] = p[i];
	    r[roff_[0]*i+2] = p[poff_[0]+i];
	    r[roff_[0]*i+0] = ((R)1.0) - (r[roff_[0]*i+1]+r[roff_[0]*i+2]);
	  } }
    }
  else if (degree[0]==0)
    {
      { I i;
	for (i=0;i<n[0];++i)
	  {
	    r[roff_[0]*i+0] = ((R)1.0);
	  } }
    }
  else if (degree[0]==2)
    {
      { I i;
	for (i=0;i<n[0];++i)
	  {
	    const R x   	= p[i];
	    const R y   	= p[poff_[0]+i];
	    const R lam 	= ((R)1.0)-(x+y);
	    r[roff_[0]*i+0] 	= lam*(((R)2.0)*lam-((R)1.0));
	    r[roff_[0]*i+1] 	= x*(((R)2.0)*x-((R)1.0));	
	    r[roff_[0]*i+2] 	= y*(((R)2.0)*y-((R)1.0));	
	    r[roff_[0]*i+3] 	= ((R)4.0)*x*lam;
	    r[roff_[0]*i+4] 	= ((R)4.0)*x*y;	  
	    r[roff_[0]*i+5] 	= ((R)4.0)*y*lam;
	  } }
    }
#if 0
  else if (degree[0]==3)
    {
      { I i;
	for (i=0;i<n[0];++i)
	  {
	    const R x   = p[i];
	    const R y   = p[poff_[0]+i];
	    const R lam = ((R)1.0)-(x+y);
	    r[roff_[0]*i+0]  = ((R)0.5)*(lam)*((-((R)1.0))+((R)3.0)*lam)*((-((R)2.0))+((R)3.0)*lam);
	    r[roff_[0]*i+1]  = ((R)0.5)*x*((-((R)1.0))+((R)3.0)*x)*((-((R)2.0))+((R)3.0)*x);
	    r[roff_[0]*i+2]  = ((R)0.5)*y*((-((R)1.0))+((R)3.0)*y)*((-((R)2.0))+((R)3.0)*y);
	    r[roff_[0]*i+3]  = ((R)0.5)*((R)9.0)*lam*x*((-((R)1.0))+((R)3.0)*lam);
	    r[roff_[0]*i+4]  = ((R)0.5)*((R)9.0)*lam*x*(((R)3.0)*x-((R)1.0));
	    r[roff_[0]*i+5]  = ((R)0.5)*((R)9.0)*x*y*(((R)3.0)*x-((R)1.0));
	    r[roff_[0]*i+6]  = ((R)0.5)*((R)9.0)*x*y*(((R)3.0)*y-((R)1.0));
	    r[roff_[0]*i+7]  = ((R)0.5)*((R)9.0)*lam*y*(((R)3.0)*y-((R)1.0));
	    r[roff_[0]*i+8]  = ((R)0.5)*((R)9.0)*lam*y*(((R)3.0)*lam-((R)1.0));
	    r[roff_[0]*i+9]  = ((R)0.5)*((R)54.0)*lam*x*y; 
	  } }
    }
#endif
  else
    {
      I info,nshapes = ((degree[0]+1)*(degree[0]+2))/2;
      I perm[nshapes];  
      mkS_lagrange_localspl_tria(degree,&rwork[0],&nshapes);
      mkS_canonic_tria(degree,&nshapes,&rwork[2*nshapes],&nshapes,&rwork[0],&nshapes,NULL,NULL,err_);
      mkS_canonic_tria(degree,n,r,roff_,p,poff_,NULL,NULL,err_);
      dgesv(&nshapes,n,&rwork[2*nshapes],&nshapes,perm,r,roff_,&info);  
    }
    
#if 0
  printf("degree= "ifmt"\n",degree[0]);
  printf("nshapes= "ifmt"\n",nshapes);
  printf("n= "ifmt"\n",n[0]);
  printf("rwork_n= "ifmt" , roff "ifmt"\n",rwork_n[0],roff_[0]);
  I i;
  for (i=0;i<nshapes*nshapes;++i) 
    printf("r=%e\n",rwork[2*nshapes+i]);
  printf("r=%e\n",r[1]);
  printf("r=%e\n",r[2]);
  printf("r=%e\n",r[3]);
  printf("r=%e\n",r[4]);
  printf("r=%e\n",r[5]);
  char A1[2] = {'A','\0'};
  char A2[2] = {'A','\0'};
  R AS[64],AU[64*64],AV[64*64],rwork2[2048];I rwork2_n=2048;
  dgesvd(A1,A2,&nshapes,&nshapes,&rwork[2*nshapes],&nshapes,AS,AU,&nshapes,AV,&nshapes,rwork2,&rwork2_n,&info);
  { I i,j;for (i=0;i<nshapes;++i){ R y =AS[i];if (nsFABS(AS[i])>((R)1.0e-15)) { y=((R)1.0)/AS[i]; for (j=0;j<nshapes;++j) {AV[j*nshapes+i] *=y; }   }else {for (j=0;j<nshapes;++j) {AV[j*nshapes+i] =(R)0.0;} } } }
  Blas_dgemm(__vmps_blas_transT,__vmps_blas_transN,&nshapes,n,&nshapes,&__vmps_blas_regal1,AU,&nshapes,r,roff_,&__vmps_blas_regal0,rwork2,&nshapes);
  Blas_dgemm(__vmps_blas_transT,__vmps_blas_transN,&nshapes,n,&nshapes,&__vmps_blas_regal1,AV,&nshapes,rwork2,&nshapes,&__vmps_blas_regal0,r,roff_);
#endif

}



