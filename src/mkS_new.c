#include "mkSYS.h"
#include "mkS.h"

mkS mkS_new(cst_emk_s kind_)
{
  emk_err err = __emk_err_no;
  mkS s = calloc(1,sizeof(mkS_st));
  mkS_def(s,
	  kind_,
	  &err);
  return s;
}

