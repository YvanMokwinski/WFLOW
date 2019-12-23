#include "nsGLOBAL.h"

pCsfInfoReadOnly 	GlobalReadOnly_get_CsfInfo(pGlobalReadOnly const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->m_csfInfo;
}
