#ifndef __header_NonLinearSolverResults_hpp__
#define __header_NonLinearSolverResults_hpp__

struct NonLinearSolverResults
{
public:
  
  R 	nrmr	[32];
  R 	nrmc	[32];
  R 	nrmh	[32];
  R 	nrmcs	[32];
  R 	nrmcu	[32];
  R 	nrmcp	[32];
  R 	nrmrs	[32];
  R 	nrmru	[32];
  R 	nrmrp	[32];
  
  NonLinearSolverResults()
  {
    Reset();
  };

  void Reset()
  {
    memset(this->nrmr,0,sizeof(R)*32);
    memset(this->nrmc,0,sizeof(R)*32);
    memset(this->nrmh,0,sizeof(R)*32);
    memset(this->nrmcs,0,sizeof(R)*32);
    memset(this->nrmcu,0,sizeof(R)*32);
    memset(this->nrmcp,0,sizeof(R)*32);
    memset(this->nrmrs,0,sizeof(R)*32);
    memset(this->nrmru,0,sizeof(R)*32);
    memset(this->nrmrp,0,sizeof(R)*32);
  };
  
};


#endif
