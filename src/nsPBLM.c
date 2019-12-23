#include "nsPBLM.h"
#include "Monitor.h"

#if 0
void use_integrator_probleme_wobbling_square	(I job_[1]);
void use_integrator_probleme_laplace		(I job_[1]);
#include "integrator_probleme_laplace.c"


void use_integrator_probleme_wobbling_square	(I job_[1]);
#include "integrator_probleme_wobbling_square.c"

void use_integrator_probleme_buckling		(I job_[1]);
#include "integrator_probleme_buckling.c"

void use_integrator_probleme_advection		(I job_[1]);
#include "integrator_probleme_advection.c"

void use_integrator_probleme_vortex		(I job_[1]);
#include "integrator_probleme_vortex.c"

void use_integrator_probleme_laplace_ellipse	(I job_[1]);
#include "integrator_probleme_laplace_ellipse.c"

void use_integrator_probleme_bench		(I job_[1]);
#include "integrator_probleme_bench.c"

void use_integrator_probleme_dambreak		(I job_[1]);
#include "integrator_probleme_dambreak.c"

void use_integrator_probleme_canal		(I job_[1]);
#include "integrator_probleme_canal.c"

void use_integrator_probleme_etirement		(const I iproc_,I job_[1]);
#include "integrator_probleme_etirement.c"

void use_integrator_probleme_rising_bubble	(const I iproc_,I job_[1]);
#include "integrator_probleme_rising_bubble.c"
#endif

void 		use_integrator_probleme_wobbling(eStateSolver job_[1]);
void 		use_integrator_probleme_canal(eStateSolver job_[1],void * );

#define table_nsPBLM_n 2
static const nsPBLM_ST table_nsPBLM[table_nsPBLM_n]
= 
  {
#if 0
    {"advection",use_integrator_probleme_advection},
    {"bench",use_integrator_probleme_bench},
    {"canal",use_integrator_probleme_canal},
    {"dambreak",use_integrator_probleme_dambreak},
    {"etirement",use_integrator_probleme_etirement},
    {"laplace",use_integrator_probleme_laplace},
    {"laplace_ellipse",use_integrator_probleme_laplace_ellipse},
    {"vortex",use_integrator_probleme_vortex},
#endif
    {"wobbling",NULL},
    {"canal",use_integrator_probleme_canal}
#if 0
    {"wobbling_square",use_integrator_probleme_wobbling_square},
    {"buckling",use_integrator_probleme_buckling},
    {"rising_bubble",use_integrator_probleme_rising_bubble},
#endif
  };


void nsPBLM_run(const char * 	name_,
		eStateSolver 	job_[1],
		void * 		usrptr_,
		STR 		errmsg_,
		Err*		err_)
{
  err_[0] = __eErr_no;
  
  { I i;
    for (i=0;i<table_nsPBLM_n;++i)
      {
	if (NOT strcmp(name_,table_nsPBLM[i].name))
	  {
	    table_nsPBLM[i].use(job_,usrptr_);
	    break;
	  }
      }

    
    if (i>=table_nsPBLM_n)
      {
	int n = sprintf(errmsg_,"'%s' is not a pblm name,\nlist of problems is : ",name_);	
	{ I j;
	  for (j=0;j<table_nsPBLM_n;++j)
	    {
	      n+=sprintf(&errmsg_[n],"\t%s\n",table_nsPBLM[j].name);	      
	    } }
	err_[0] = __eErr_user;
	return;
      }
    
  }
}

