#include "nsGLOBAL.h"


pVariables			Global_get_Variables		(pGlobal const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return &self_->m_variables;
}
