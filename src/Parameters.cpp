#include "Parameters.hpp"

typedef void * pParameters;
typedef const void * pParametersReadOnly;
extern "C"
{
  L 			Parameters_getl			(const void* 	self_,
							 const ens_linfo linfo_)
  {
    const Parameters * p = (const Parameters*)self_;
    return p->GetInfoLogical((InfoLogical::EnumType)linfo_) ? __emnsYES: __emnsNO;
  };

  R 			Parameters_getr			(const void* 	self_,
							 const ens_rinfo rinfo_)
  {
    const Parameters * p = (const Parameters*)self_;
    return p->GetInfoReal((InfoReal::EnumType)rinfo_);

  };
  
  I 			Parameters_geti			(const void* 	self_,
							 const ens_iinfo iinfo_)
  {
    const Parameters * p = (const Parameters*)self_;
    return p->GetInfoInteger((InfoInteger::EnumType)iinfo_);

  };

  eTransientMethod 	Parameters_get_eTransientScheme	(const void*self_,
							 eKindEquation	kindEquation_)
  {
    const Parameters * p = (const Parameters*)self_;
    return (eTransientMethod)p->GetKindTransientMethod( (KindEquation::EnumType) kindEquation_);
  };

  
}




#if 0
eLinearSolver 	Parameters_get_eLinearSolver	(const void*self_,
						 eKindEquation	kindEquation_)
{

};

extern "C"
{

eLinearSolver 	ParametersReadOnly_get_eLinearSolver	(pParametersReadOnly 	const 	self_,
							 eKindEquation 			kindEquation_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->method_linearsolver[kindEquation_];
}

L ParametersReadOnly_getl(pParametersReadOnly const self_,const ens_linfo linfo_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->linfo[linfo_];
}

void Parameters_free(pParameters const params_)
{
#ifndef NDEBUG
  DebugVerif(params_);
#endif
  if (params_)
    {
      memset(params_,0,sizeof(Parameters));
    }
}

Err Parameters_def(pParameters 	const 	params_,
		   pCmdline 	const 	cmdline_,
		   const L  		have_configfile_)
{
#ifndef NDEBUG
  DebugVerif(params_);
  DebugVerif(cmdline_);
#endif
  ens_linfo_from_Cmdline(params_->linfo,have_configfile_,cmdline_);
  ens_iinfo_from_Cmdline(params_->iinfo,have_configfile_,cmdline_);
  ens_rinfo_from_Cmdline(params_->rinfo,have_configfile_,cmdline_);
  ens_sinfo_from_Cmdline(params_->sinfo,have_configfile_,cmdline_);
  return __eErr_no;
}

};
#endif
