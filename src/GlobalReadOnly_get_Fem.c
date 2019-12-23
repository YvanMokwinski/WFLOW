#include "nsGLOBAL.h"


pFemReadOnly	GlobalReadOnly_get_Fem	(pGlobalReadOnly const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->m_fem;
}
