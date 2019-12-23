#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Cmdline.h"
#include "Monitor.h"

void Cmdline_disp(cst_pCmdline 	self_)
{
  { int i; 
    for (i=0;i<self_->argc;++i) 
      { 
	Monitor_msg(0,"[%d] %s\n",i,self_->argv[i]);
      } }
}


void Cmdline_free	(pCmdline 		cmdline_)
{
  if (cmdline_)
    {
      memset(cmdline_,0,sizeof(Cmdline));
    }
}

void Cmdline_def	(pCmdline 		cmdline_,
			 const int		argc_,
			 cst_pS * 		argv_)
{
 
  if (argc_<__mCmdline_MAXARG)
    {
      memset(cmdline_,0,sizeof(Cmdline));
      cmdline_->ref_argc 	= argc_;
      cmdline_->ref_argv 	= argv_;
      cmdline_->argc 	= argc_;
      { int i; 
	for (i=0;i<argc_;++i) 
	  { 
	    strncpy( cmdline_->argv[i],argv_[i],(size_t)__mSTR_maxlen); 
	  } }
    }
  else
    {
      /* ns_errmsg("Cmdline_def:__mCmdline_MAXARG=%d is reached",argc_);*/
    }
}


static int Cmdline_search(int 		argc_,
			  STR* 		argv_,
			  cst_pS 	opt_,
			  const I nchoice_,
			  cst_pS* 	opt_choice_)
{
  int j,i = (int)0;
  for (i=1;i<argc_;i++)
    if (!strcmp(argv_[i],opt_))
      break;
  if (i<argc_)
    {
      if ( (nchoice_>0) && (i+1<argc_) )
	{
	  ++i;
	  for (j=0;j<nchoice_;++j)
	    if (! strcmp(argv_[i],opt_choice_[j]))
	      break;
	  if (j<nchoice_)
	    return i-1;
	  else 
	    return (int)-1;
	}      
      else
	return i;
    }
  else
    return (int)-1;
}

L Cmdline_get_logical(pCmdline 	cmdline_,
		      const STR	opt_)
{  
  int k = Cmdline_search(cmdline_->argc,
			 cmdline_->argv,
			 opt_,
			 (I)0,
			 NULL);
  
  if ( (k>=0) && (k<cmdline_->argc) )
    {
      if (k+1<cmdline_->argc)
	{
	  { int i;
	    for (i=0;i<cmdline_->argc-k;++i)
	      {
		strncpy(cmdline_->argv[k+i],cmdline_->argv[k+1+i],(size_t)__mSTR_maxlen);
	      } }
	}
      cmdline_->argc-=1;
      return __emnsYES;
    }
  return __emnsNO;
}


L Cmdline_get_string(pCmdline 	cmdline_,
			       const STR 		opt_,
			       STR 			t_)
{  
  int j,k = Cmdline_search(cmdline_->argc,
			     cmdline_->argv,
			     opt_,
			     (I)0,
			     NULL);
  if ( (k>=0) && (++k<cmdline_->argc) )
    {
      j=k-1;
      strncpy(t_,cmdline_->argv[k],(size_t)__mSTR_maxlen);
      {
	int i;
	for (i=0;i<cmdline_->argc-k;++i)
	  {
	    strncpy(cmdline_->argv[j+i],cmdline_->argv[k+1+i],(size_t)__mSTR_maxlen);
	  }
      }
      cmdline_->argc-=2;
      return __emnsYES;
    }
  return __emnsNO;
}

L  Cmdline_get_integer(pCmdline 	cmdline_,
			       const STR	opt_,
			       pI 		x_)
{  
  int j,k = Cmdline_search(cmdline_->argc,
			     cmdline_->argv,
			     opt_,
			     (I)0,
			     NULL);
  if ( (k>=0) && (++k<cmdline_->argc) )
    {
      j=k-1;
      sscanf(cmdline_->argv[k],""ifmt"",x_);
      {	int i;
	for (i=0;i<cmdline_->argc-k;++i)
	  {
	    strncpy(cmdline_->argv[j+i],cmdline_->argv[k+1+i],(size_t)__mSTR_maxlen);
	  } }
      cmdline_->argc-=2;
      return __emnsYES;
    }
  return __emnsNO;
}


L  Cmdline_get_real(pCmdline 	cmdline_,
		    const STR 		opt_,
		    pR 			x_)
{  
  int j,k = Cmdline_search(cmdline_->argc,
			     cmdline_->argv,
			     opt_,
			     (I)0,
			     NULL);
  if ( (k>=0) && (++k<cmdline_->argc) )
    {
      j=k-1;
      sscanf(cmdline_->argv[k],"%"rformat"",x_);
      {	int i;
	for (i=0;i<cmdline_->argc-k;++i)
	  {
	    strncpy(cmdline_->argv[j+i],cmdline_->argv[k+1+i],(size_t)__mSTR_maxlen);
	  } }
      cmdline_->argc-=2;
      return __emnsYES;
    }
  return __emnsNO;
}

Err Cmdline_check_invalid(cst_pCmdline cmdline_)
{
  { int i,j=0;
    for (i=2;i<cmdline_->argc;++i)
      {
	if (cmdline_->argv[i][0]=='-')
	  {
	    j=1;
	    /*	    ns_errmsg("main:unknown option '%s'",cmdline_->argv[i]);*/
	  }
      } 
    if (j>0)
      {
	return __eErr_user;
      } }  
  return __eErr_no;
}


L Cmdline_isempty	(cst_pCmdline cmdline_)
{
  return (cmdline_->argc==1)?__emnsYES:__emnsNO;
}



cst_pS	Cmdline_get_arg	(pCmdline 	cmdline_,
			 const int 	i_)
{
  return (i_<cmdline_->argc) ? (cmdline_->argv[i_]) : NULL;
}



