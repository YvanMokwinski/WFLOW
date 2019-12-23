#include "nsGLOBAL.h"


pCsfInfo 		Global_get_CsfInfo	(pGlobal const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->m_csfInfo;
}
