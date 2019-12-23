#include "Err.h"
#include <stdio.h>
#include <stdlib.h>

cst_pS eErr_get_string(const Err err_)
{
  cst_pS str = 0;
  switch(err_)
    {

    case __eErr_no:
      {
	str = "NO ERROR";
	break;
      }

    case __eErr_user:
      {
	str = "USER ERROR";
	break;
      }

    case __eErr_file:
      {
	str = "FILE ERROR";
	break;
      }

    case __eErr_switch:
      {
	str = "SWITCH ERROR";
	break;
      }

    case __eErr_blank:
      {
	str = "BLANK ERROR";
	break;
      }

    case __eErr_iwork:
      {
	str = "IWORK ERROR";
	break;
      }

    case __eErr_intern:
      {
	str = "INTERN ERROR";
	break;
      }

    case __eErr_memory:
      {
	str = "MEMORY ERROR";
	break;
      }

    case __eErr_ALL:
      {
	str = "ALL ERRORS ?";
	break;

      }
    }
  return str;
}

void   eErr_check	(const Err 	err_,
			 cst_pS 	file_,
			 const I 	line_)
{
  if (err_)
    {
      cst_pS str = eErr_get_string(err_);
      fprintf(stderr,"ERROR (%s) from file '%s', line "ifmt"\n",str,file_,line_);
      exit(1);
    }
}
