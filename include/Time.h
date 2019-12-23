#ifndef __header_Time_h__
#define __header_Time_h__

#include "Type.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct
  {
    R 		t;
    R 		ti;
    R 		tii;
    R 		t0;
    R 		tf;
    R 		time_interval_t0;
    R 		time_interval_tf;
    R 		time_interval_tm2;
    R 		time_interval_tm1;
    I 		itime;	
    I 		itime_interval;
    I 		itime_interval_tm1;
    I 		itime_interval_tm2;  
  } Time,* RESTRICT pTime;

  typedef const Time * RESTRICT 	pTimeReadOnly;

  I 	TimeReadOnly_get_itime			(pTimeReadOnly 	const 	self_);
  R 	TimeReadOnly_get_t			(pTimeReadOnly 	const 	self_);
  R 	TimeReadOnly_get_ti			(pTimeReadOnly 	const 	self_);
  R 	TimeReadOnly_get_tii			(pTimeReadOnly 	const 	self_);
  void 	TimeReadOnly_get_allt			(pTimeReadOnly 	const 	self_,
						 pR 			allt_);
  R 	TimeReadOnly_get_dt			(pTimeReadOnly 	const 	self_);
  R 	TimeReadOnly_get_dti			(pTimeReadOnly 	const 	self_);

  void 	Time_set_allt			(pTime 		const 	self_,
					 cst_pR			allt_);

  void 	Time_restart_time_interval	(pTime		const	self_);
  void 	Time_next_time_interval		(pTime		const	self_,
					 const R 		a_);
  void 	Time_next_interval		(pTime		const	self_);
  void 	Time_next_time			(pTime		const	self_,
					 const R 		a_);
  void 	Time_init			(pTime		const	self_);


#ifdef __cplusplus
}
#endif

#endif
