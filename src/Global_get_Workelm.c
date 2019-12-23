#include "nsGLOBAL.h"

pWorkelm Global_get_Workelm(pGlobal 		const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->m_workelm;
}
