#include "nsGLOBAL.h"

pTime			Global_get_Time		(pGlobal const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->m_timeInfo;
}
