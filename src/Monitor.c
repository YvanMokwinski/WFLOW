#include "Monitor.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "Config.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/**
   \brief Type array of character
*/
typedef char STR[256];
/**
   \brief Declaration of static variables to handle files, here's mode for monitoring.
*/
static unsigned char 	modes		[MNS_MAX_THREAD] = {MonitorMode_STD};
/**
   \brief Declaration of static variables to handle files, here's string for program names.
*/
static STR 		prognames	[MNS_MAX_THREAD] = { {0} };
/**
   \brief Declaration of static variables to handle files, here's pointer for files.
*/
static FILE* 		logfiles	[MNS_MAX_THREAD] = {0};


#ifndef NDEBUG
/* is undefined with 'undef' macro at the end of this file */
#define DebugVerif(_cond) if (!(_cond)) { fprintf(stderr,"condition failed line %d file %s\n",__LINE__,__FILE__);   }
#endif


void Monitor_createDirectory(const char * 	name_,
			     ...)
{
  struct stat st = {0};
#if 0
  STR namebis;    sprintf(namebis,"%s.%d.%s",prognames[ithread_],ithread_,name);  
#endif

  STR name;  
  { va_list args;
    va_start (args,name_);
    vsprintf(name,name_,args);
    va_end(args); }


  if (stat(name, &st) == -1) 
    {
      mkdir(name, 0700);
    }
}


/**/

void Monitor_shutdown()
{
  { int i;
    for (i=0;i<MNS_MAX_THREAD;++i)
      {
	if (logfiles[i])
	  fclose(logfiles[i]);
      } }
}

/**/

void Monitor_def(const int 		ithread_,
		 const char*		progname_,
		 const unsigned char 	mode_)
{
#ifndef NDEBUG
  DebugVerif(ithread_<MNS_MAX_THREAD);
  DebugVerif(mode_>0);
#endif


  sprintf	(prognames[ithread_],
		 "%s",
		 progname_);
  modes[ithread_]	= mode_;
  unsigned char mode = modes[ithread_];
  mode = mode >> 1;
  if (mode%2>0)
    {
      STR a;
      const char * progname_without_path 	= NULL;
      const char * tmp 				= NULL;
      for (tmp=prognames[ithread_];tmp[0]!='\0';++tmp)
	if (*tmp=='/')
	  progname_without_path = tmp;
      progname_without_path = (progname_without_path)?progname_without_path+1:prognames[ithread_];
      sprintf(a,"%s.thread%d.log",progname_without_path,ithread_);
      logfiles[ithread_] = fopen(a,"w");
    }



}

/**/

void Monitor_msg(const int 		ithread_,
		 const char* 		msg_,
		 ...) 
{
#ifndef NDEBUG
  DebugVerif(msg_);
  DebugVerif(ithread_<MNS_MAX_THREAD);
#endif
  STR  		ctmp;
  unsigned char mode = modes[ithread_];

  { va_list args;
    va_start (args,msg_);
    vsprintf(ctmp,msg_,args);
    va_end(args); }

  if (mode%2>0)
    {
      /*MonitorMode_STD*/
      fflush(stderr);
      fprintf(stderr,"//%s.%d:MSG:%s\n",prognames[ithread_],ithread_,ctmp);  
    }

  mode = mode >> 1;
  if (mode%2>0)
    {
      /*MonitorMode_LOG*/
      fprintf(logfiles[ithread_],"//%s.%d:MSG:%s\n",prognames[ithread_],ithread_,ctmp);  
    }
}

/**/

void Monitor_errmsg(const int ithread_,
		    const char* msg_,
		    ...) 
{
#ifndef NDEBUG
  DebugVerif(msg_);
  DebugVerif(ithread_<MNS_MAX_THREAD);
#endif
  STR  		ctmp;
  unsigned char mode = modes[ithread_];

  { va_list args;
    va_start (args,msg_);
    vsprintf(ctmp,msg_,args);
    va_end(args); }

  if (mode%2>0)
    {
      /*MonitorMode_STD*/
      fflush(stderr);
      fprintf(stderr,"//%s.%d:ERROR:%s\n",prognames[ithread_],ithread_,ctmp);  
    }

  mode = mode >> 1;
  if (mode%2>0)
    {
      /*MonitorMode_LOG*/
      fprintf(logfiles[ithread_],"//%s.%d:ERROR:%s\n",prognames[ithread_],ithread_,ctmp);  
    }

  mode = mode >> 1;
  if (mode%2>0)
    {
      /*MonitorMode_EXITIFERROR*/
      exit(1);
    }
}

/**/

void Monitor_warn(const int ithread_,
		  const char* msg_,
		  ...) 
{
#ifndef NDEBUG
  DebugVerif(msg_);
  DebugVerif(ithread_<MNS_MAX_THREAD);
#endif
  STR  		ctmp;
  unsigned char mode = modes[ithread_];

  { va_list args;
    va_start (args,msg_);
    vsprintf(ctmp,msg_,args);
    va_end(args); }


  if (mode%2>0)
    {
      /*MonitorMode_STD*/
      fflush(stderr);
      fprintf(stderr,"//%s.%d:WARN:%s\n",prognames[ithread_],ithread_,ctmp);  
    }

  mode = mode >> 1;
  if (mode%2>0)
    {
      /*MonitorMode_LOG*/
      fprintf(logfiles[ithread_],"//%s.%d:WARN:%s\n",prognames[ithread_],ithread_,ctmp);  
    }

}

/**/

#ifndef NDEBUG
#undef DebugVerif
#endif
