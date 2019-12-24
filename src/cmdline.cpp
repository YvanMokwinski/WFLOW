
#include "cmdline.hpp"
#include <string.h>  
#include <stdio.h>  

//!
//!   @brief Class to define a command line.
//!
int cmdline::search(int 		argc_,
		    char argv_[cmdline::MAXARG][cmdline::MAXLEN] ,
		    const char* 	opt_,
		    const int 	nchoice_,
		    const char** opt_choice_)
{
  int j,i = (int)0;
  for (i=1;i<argc_;i++)
    {
      if (!strcmp(argv_[i],opt_))
	{
	  break;
	}
    }
  if (i<argc_)
    {
      if ( (nchoice_>0) && (i+1<argc_) )
	{
	  ++i;
	  for (j=0;j<nchoice_;++j)
	    {
	      if (! strcmp(argv_[i],opt_choice_[j]))
		{
		  break;
		}
	    }
	  if (j<nchoice_)
	    {
	      return i-1;
	    }
	  else 
	    {
	      return (int)-1;
	    }
	}      
      else
	{
	  return i;
	}
    }
  else
    {
      return (int)-1;
    }
};

void cmdline::disp() const
{

  for (int i=0;i<this->argc;++i) 
    { 
      fprintf(stdout,"[%d] %s\n",i,this->argv[i]);
    } 
};


cmdline::~cmdline()
{
};

cmdline::cmdline(int		argc_,
		 char ** 		argv_)
{  
  if (argc_<cmdline::MAXARG)
    {
      this->ref_argc 	= argc_;
      this->ref_argv 	= argv_;
      this->argc 	= argc_;
      for (int i=0;i<argc_;++i) 
	{ 
	  strncpy( this->argv[i],argv_[i],(size_t)cmdline::MAXLEN); 
	} 
    }
  else
    {
      // ns_errmsg("cmdline::def:__mcmdline::MAXARG=%d is reached",argc_);
    }
};



bool cmdline::get_logical(const char*	opt_)
{  
  int k = cmdline::search(this->argc,
			  this->argv,
			  opt_);
  
  if ( (k>=0) && (k<this->argc) )
    {
      if (k+1<this->argc)
	{
	  { int i;
	    for (i=0;i<this->argc-k;++i)
	      {
		strncpy(this->argv[k+i],this->argv[k+1+i],(size_t)cmdline::MAXLEN);
	      } }
	}
      this->argc-=1;
      return true;
    }
  return false;
};


bool cmdline::get_string(const char* 		opt_,
			 char* 			t_)
{  
  int j,k = cmdline::search(this->argc,
			    this->argv,
			    opt_);
  if ( (k>=0) && (++k<this->argc) )
    {
      j=k-1;
      strncpy(t_,this->argv[k],(size_t)cmdline::MAXLEN);
      {
	int i;
	for (i=0;i<this->argc-k;++i)
	  {
	    strncpy(this->argv[j+i],this->argv[k+1+i],(size_t)cmdline::MAXLEN);
	  }
      }
      this->argc-=2;
      return true;
    }
  return false;
};



template<typename int_t>
bool  cmdline::get_integer(const char*	opt_,
			   int_t* 	x_)
{  
  int j,k = cmdline::search(this->argc,
			     this->argv,
			    opt_);
  if ( (k>=0) && (++k<this->argc) )
    {
      j=k-1;
      sscanf(this->argv[k],"%d",x_);
      {	int i;
	for (i=0;i<this->argc-k;++i)
	  {
	    strncpy(this->argv[j+i],this->argv[k+1+i],(size_t)cmdline::MAXLEN);
	  } }
      this->argc-=2;
      return true;
    }
  return false;
};

template<>
bool  cmdline::get_integer<long long int>(const char*	opt_,
				      long long int* 	x_)
{  
  int j,k = cmdline::search(this->argc,
			     this->argv,
			    opt_);
  if ( (k>=0) && (++k<this->argc) )
    {
      j=k-1;
      sscanf(this->argv[k],"%Ld",x_);
      {	int i;
	for (i=0;i<this->argc-k;++i)
	  {
	    strncpy(this->argv[j+i],this->argv[k+1+i],(size_t)cmdline::MAXLEN);
	  } }
      this->argc-=2;
      return true;
    }
  return false;
};


template<typename real_t>
bool  cmdline::get_real(const char* 		opt_,
			real_t* 		x_)
{  
  int j,k = cmdline::search(this->argc,
			    this->argv,
			    opt_);
  if ( (k>=0) && (++k<this->argc) )
    {
      j=k-1;
      sscanf(this->argv[k],"%le",x_);
      {	int i;
	for (i=0;i<this->argc-k;++i)
	  {
	    strncpy(this->argv[j+i],this->argv[k+1+i],size_t(cmdline::MAXLEN));
	  } }
      this->argc-=2;
      return true;
    }
  return false;
};


template<>
bool  cmdline::get_real<double>(const char* 		opt_,
				double* 		x_)
{  
  int j,k = cmdline::search(this->argc,
			    this->argv,
			    opt_);
  if ( (k>=0) && (++k<this->argc) )
    {
      j=k-1;
      sscanf(this->argv[k],"%le",x_);
      {	int i;
	for (i=0;i<this->argc-k;++i)
	  {
	    strncpy(this->argv[j+i],this->argv[k+1+i],size_t(cmdline::MAXLEN));
	  } }
      this->argc-=2;
      return true;
    }
  return false;
};


bool cmdline::check_invalid() const
{
  int i,j=0;
  for (i=2;i<this->argc;++i)
    {
      if (this->argv[i][0]=='-')
	{
	  j=1;
	}
    }
  return j>0; 
};


bool cmdline::isempty() const
{
  return (1 == this->argc);
};

const char*	cmdline::get_arg(const int i_) const 
{
  return (i_<this->argc) ? (this->argv[i_]) : nullptr;
};
