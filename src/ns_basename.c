#include "ns_config.h"
#include "Type.h"

void ns_basename(const char * 	t_,
		 char		o_[512])
{
  I j=0;
  while((j<512)AND(t_[j]!='\0'))++j;
  I n = ((j<512)?j:0);j=0;{I i;for(i=0;i<n;++i) if(t_[i]=='.')j=i;}const char*pp=(!j)?NULL:(&t_[j+1]);if(!pp)j=n;{I i;for(i=0;i<j;++i){o_[i]=t_[i];}}o_[j]='\0';
}

