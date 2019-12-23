#include "Fem.h"

cst_pR  FemReadOnly_get_dens_u_u	(pFemReadOnly	const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->dens_u_u;
}

pR  	Fem_get_dens_u_u		(pFem	const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->dens_u_u;
}



cst_pR  FemReadOnly_getJ_UU		(pFemReadOnly	const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->J_UU;
}

pR  	Fem_getJ_UU			(pFem		const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->J_UU;
}


cst_pR  FemReadOnly_getJ_VV		(pFemReadOnly	const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->J_VV;
}

pR  	Fem_getJ_VV			(pFem		const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->J_VV;
}

cst_pR  FemReadOnly_getJ_UV		(pFemReadOnly	const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->J_UV;
}

pR  	Fem_getJ_UV			(pFem		const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->J_UV;
}

cst_pR  FemReadOnly_getJ_VU		(pFemReadOnly	const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->J_VU;
}

pR  	Fem_getJ_VU			(pFem		const 	self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->J_VU;
}
