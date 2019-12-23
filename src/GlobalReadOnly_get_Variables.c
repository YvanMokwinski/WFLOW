#include "nsGLOBAL.h"


pVariablesReadOnly			GlobalReadOnly_get_Variables		(pGlobalReadOnly const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->m_variables;
}
