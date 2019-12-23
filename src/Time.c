#include "Time.h"

R TimeReadOnly_get_dt(pTimeReadOnly const self_) 
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->t-self_->ti;
}

R TimeReadOnly_get_dti(pTimeReadOnly const self_) 
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->ti-self_->tii;
}

I TimeReadOnly_get_itime(pTimeReadOnly const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->itime;
}

I TimeReadOnly_get_itime_interval(pTimeReadOnly const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->itime_interval;
}

R TimeReadOnly_get_t(pTimeReadOnly const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->t;
}

R TimeReadOnly_get_ti(pTimeReadOnly const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->ti;
}

R TimeReadOnly_get_tii(pTimeReadOnly const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  return self_->tii;
}


void TimeReadOnly_get_allt(pTimeReadOnly const self_,
		   pR allt_)
{
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(allt_);
#endif
  allt_[0] = self_->t;
  allt_[1] = self_->ti;
  allt_[2] = self_->tii;
}


void Time_set_allt(pTime const self_,
		   cst_pR allt_)
{
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(allt_);
#endif
  self_->t 	= allt_[0];
  self_->ti 	= allt_[1];
  self_->tii 	= allt_[2];
}

void Time_restart_time_interval(pTime const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  self_->t     		= self_->time_interval_t0;
  self_->ti    		= self_->time_interval_tm1;
  self_->tii   		= self_->time_interval_tm2;
}

void Time_next_time_interval(pTime const self_,const R a_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  self_->time_interval_t0 	= self_->time_interval_tf;
  self_->time_interval_tf 	= self_->time_interval_tf+a_;
}


void Time_next_interval(pTime const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  ++self_->itime_interval;
  self_->time_interval_tm2 	= self_->tii;
  self_->time_interval_tm1 	= self_->ti;
  self_->itime_interval_tm1 	= self_->itime-1;
  self_->itime_interval_tm2 	= self_->itime-2;
}

void Time_next_time(pTime const self_,const R a_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  self_->tii   			= self_->ti;
  self_->ti    			= self_->t;
  self_->t     			= self_->ti;
}

void Time_init(pTime const self_)
{
#ifndef NDEBUG
  DebugVerif(self_);
#endif
  self_->time_interval_t0 	= self_->t0;
  self_->time_interval_tf 	= self_->t0;
  self_->time_interval_tm2 	= self_->t0;
  self_->time_interval_tm1 	= self_->t0;
  self_->itime_interval 	= 0;
  self_->itime 			= 0;
}
