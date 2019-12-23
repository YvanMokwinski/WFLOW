#include "nsGLOBAL.h"
#include "ns_config.h"
#include "libns_exact.h"
#include "Monitor.h"

void Global_precompute(pGlobal 		const 	self_,
		       ns_mesh * 	const 	mesh_,
		       STR 			errmsg_,
		       Err* 			err_)
{
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(mesh_);
#endif
  err_[0] 			= __eErr_no;
  pParameters 	parameters 	= Global_get_Parameters(self_);
  pWorkelm	workelm		= Global_get_Workelm(self_);
  const L 	color 		= ParametersReadOnly_getl(parameters,__ens_linfo_color);
  
  /*1 ###############################################*/
  
#ifndef NDEBUG
  const L verbose = ParametersReadOnly_getl(parameters,__ens_linfo_verbose);
  if (verbose)
    Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 1/5:reading mesh from file '%s'",self_->filename);
#endif

  self_->mesh_usrptr = mesh_;
  Workelm_init(workelm);

  //
  // 2
  //
#ifndef NDEBUG
  if (verbose)
    Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 2/5:initialize exact data ...");
#endif
  
  Fem_def(Global_get_Fem(self_),
	  mesh_->elm,
	  errmsg_,
	  err_);
  
  if (err_[0])
    {
      Monitor_errmsg(self_->iproc,"nsGLOBAL_precompute:ns_exact_init failed");
      return;
    }

#ifndef NDEBUG
  if (verbose)
    {
      Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 2/5:initialize exact data done");
    }
#endif

  //
  // 3
  //

#ifndef NDEBUG
  if (verbose)
    {
      Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 3/5:nsGLOBAL_initialize ...");
    }
#endif

  STR errmsg;
  Global_initialize(self_,
		    (ns_mesh*)self_->mesh_usrptr,
		    errmsg,
		    err_);
  if (err_[0])
    {
      Monitor_errmsg(self_->iproc,"nsGLOBAL_precompute:nsGLOBAL_initialize failed");
      return;
    }
  
#ifndef NDEBUG
  if (verbose)
    {
      Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 3/5:nsGLOBAL_initialize done.");
    }
#endif

  //
  // 4
  //
  
#ifndef NDEBUG
  if (verbose)
    {
      Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 4/5:nsCSF_INFO_def ...");
    }
#endif

  CsfInfo_def	(&self_->m_csfInfo,
		 color,
		 __eOperator_ih,
		 __eOperator_qh,
		 __eOperator_ih,
		 ((R)0.125),
		 __eOperator_qh,
		 ((R)0.99),
		 __eOperator_ph,
		 __eOperator_ph,
		 0.5,
		 __ens_method_capturing_linear,
		 1,
		 __eSmoothingMethod_LQ_l,
		 1,
		 (eHeaviside)parameters->iinfo[__ens_iinfo_heaviside],
		 (eDirac)parameters->iinfo[__ens_iinfo_dirac],
		 &self_->rv[__ens_rv_eps],
		 err_);  
  
#ifndef NDEBUG
  if (verbose)
    {
      Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 4/5:CsfInfo_def done.");
      Monitor_msg(self_->iproc,"nsGLOBAL_precompute:step 5/5:completed.");
    }
#endif
  
}
