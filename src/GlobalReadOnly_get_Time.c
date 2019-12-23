#include "nsGLOBAL.h"


pTimeReadOnly		GlobalReadOnly_get_Time	(pGlobalReadOnly const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->m_timeInfo;
}
