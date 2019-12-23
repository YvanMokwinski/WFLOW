#include "mkSYS.h"
#include "mkS.h"

void mok_s_errmsg(const char * msg,...)
{
  char 	roger[512];
  va_list args;
  va_start (args,msg);
  fflush(stdout);
  fprintf(stderr,"\n**** ERROR mok_s.%d ****\n**** msg : ",getpid());
  vsprintf(roger,msg,args);
  fprintf(stderr,"%s\n",roger);
  va_end(args);
}
