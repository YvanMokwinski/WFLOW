#ifndef __header_Parameters_h__
#define __header_Parameters_h__

// #include "Cmdline.h"
#include "ns_enum.h"

#ifdef __cplusplus
extern "C"
{
#endif


  typedef void * pParameters;
  typedef const void * pParametersReadOnly;
  
  L 			Parameters_getl			(const void* 	self_,
							 const ens_linfo linfo_);
  R 			Parameters_getr			(const void* 	self_,
							 const ens_rinfo rinfo_);
  I 			Parameters_geti			(const void* 	self_,
							 const ens_iinfo iinfo_);


  eTransientMethod 	Parameters_get_eTransientScheme	(const void*self_,
							 eKindEquation	kindEquation_);

#if 0
  
  typedef void * pParameters;
  typedef const void * pParametersReadOnly;
  
  L 			ParametersReadOnly_getl			(pParametersReadOnly 	const 	self_,
								 const ens_linfo 		linfo_);
  
  eLinearSolver 	ParametersReadOnly_get_eLinearSolver	(pParametersReadOnly 	const 	self_,
								 eKindEquation 			kindEquation_);
#if 0  
  Err 			Parameters_def				(pParameters 		const 	self_,
								 pCmdline		const 	cmdline_,
								 const L			have_configfile);
#endif  
  void 			Parameters_free				(pParameters  		const 	self_);
#endif
  
#ifdef __cplusplus
}
#endif

#endif
