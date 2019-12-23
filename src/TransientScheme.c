
#include "TransientScheme.h"
#include "Blas.h"


Err TransientScheme_def(pTransientScheme const	scheme_,
			cst_eTransientMethod 	method_)
{
  Err err = __eErr_no;
  memset(scheme_,0,sizeof(TransientScheme));
  scheme_->method  = method_;
  switch(method_)
    {
    case __eTransientMethod_NO:
      {
	scheme_->nstep   = 0;
	break;
      }
    case __eTransientMethod_EULER:
      {
	scheme_->nstep   = 1;
	break;
      }
    case __eTransientMethod_IMR:
      {
	scheme_->nstep   = 1;
	break;
      }
    case __eTransientMethod_GEAREULER:
      {
	scheme_->nstep   = 2;
	break;
      }
    case __eTransientMethod_GEARIMR:
      {
	scheme_->nstep   = 2;
	break;
      }
    case __eTransientMethod_TRAPEZE:
      {
	scheme_->nstep   = 1;
	break;
      }
    case __eTransientMethod_ERROR:
    case __eTransientMethod_ALL:
      {
	err = __eErr_switch;
	break;
      }
    }
  if (err)
    {
      memset(scheme_,0,sizeof(TransientScheme));
      return err;
    }
  return err;
}


Err TransientScheme_precompute(pTransientScheme 	const	scheme_,
			       const I				itime_,
			       const R 				dt_,
			       const R 				idt_,
			       const R 				dti_,
			       const R 				ratio_dti_,
			       const R 				iratio_dti_)
{
  Err err = __eErr_no;
  switch(scheme_->method)
    {
    case __eTransientMethod_NO:
      {
	scheme_->scheme[0] = ((R)0.0);
	scheme_->nstep = 0;
	break;
      }
    case __eTransientMethod_EULER:
      {
	scheme_->scheme	[0] = idt_;
	scheme_->scheme	[1] = -idt_;
	scheme_->nstep = 1;
	break;
      }
    case __eTransientMethod_IMR:
      {
	scheme_->scheme	[0] = ((R)2.0)*idt_;
	scheme_->scheme	[1] = -((R)2.0)*idt_;
	scheme_->nstep = 1;
	break;
      }
    case __eTransientMethod_GEAREULER:
      {
	if (itime_>0)
	  {
	    scheme_->scheme	[0] = (((R)2.0)+ratio_dti_)/(dt_+dti_);
	    scheme_->scheme	[1] = -(((R)2.0) + ratio_dti_ + iratio_dti_)/(dt_+dti_);
	    scheme_->scheme	[2] = iratio_dti_/(dt_+dti_);
	    scheme_->nstep 	= 2;

	  }
	else
	  {
	    scheme_->scheme	[0] = idt_;
	    scheme_->scheme	[1] = -idt_;
	    scheme_->nstep 	= 1;
	  }
	break;
      }
    case __eTransientMethod_GEARIMR:
      {
	if (itime_>0)
	  {
	    scheme_->scheme	[0] = (((R)2.0)+ratio_dti_)/(dt_+dti_);
	    scheme_->scheme	[1] = -(((R)2.0) + ratio_dti_ + iratio_dti_)/(dt_+dti_);
	    scheme_->scheme	[2] = iratio_dti_/(dt_+dti_);
	    scheme_->nstep 	= 2;
	  }
	else
	  {
	    scheme_->scheme	[0] = ((R)2.0)*idt_;
	    scheme_->scheme	[1] = -((R)2.0)*idt_;
	    scheme_->nstep 	= 1;
	  }	
	break;
      }
    case __eTransientMethod_TRAPEZE:
      {
	scheme_->scheme	[0] = idt_*((R)0.5);
	scheme_->scheme	[1] = idt_*((R)0.5);
	scheme_->nstep 	= 1;
	break;
      }
    case __eTransientMethod_ERROR:
    case __eTransientMethod_ALL:
      {
	err = __eErr_switch;
	break;
      }
    }
  return err;
}

void TransientScheme_compute(pTransientScheme 	const 	scheme_,
			     cst_pI		const	fn_,
			     cst_pR		const	f_,
			     cst_pI		const	foff_,
			     pR			const	r_)
{
  static const R 	regal1		= ((R)1.0);
  static const R 	regal0		= ((R)0.0);
  static const I  	negal1		= ((I)1);
  static const char 	transN[2]	= {'N','\0'};
  const I n 				= scheme_->nstep+1;
  Blas_dgemv(transN,fn_,&n,&regal1,f_,foff_,scheme_->scheme,&negal1,&regal0,r_,&negal1);
}
