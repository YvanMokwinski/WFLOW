#include "Workelm.h"


void Workelm_init(pWorkelm const self_)
{
  self_->uelm 	= &self_->ddlelm[_ju];
  self_->velm 	= &self_->ddlelm[_jv];
  self_->pelm 	= &self_->ddlelm[_jp];

  self_->uelmi 	= &self_->ddlelmi[_ju];
  self_->velmi 	= &self_->ddlelmi[_jv];
  self_->pelmi 	= &self_->ddlelmi[_jp];

  self_->uelmii = &self_->ddlelmii[_ju];
  self_->velmii = &self_->ddlelmii[_jv];
  self_->pelmii = &self_->ddlelmii[_jp];
}
