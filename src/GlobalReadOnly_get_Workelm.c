#include "nsGLOBAL.h"

pWorkelmReadOnly GlobalReadOnly_get_Workelm(pGlobalReadOnly const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->m_workelm;
}
