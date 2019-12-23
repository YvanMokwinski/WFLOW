#include "nsGLOBAL.h"


pFem			Global_get_Fem		(pGlobal const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->m_fem;
}
