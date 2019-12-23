#ifndef __header__eStateSolver_h__
#define __header__eStateSolver_h__

#include "ns_sys.h"
#include "ConfigEnum.h"

ConfigEnum(eStateSolver,
	   __eStateSolver_dirichlet_pressure,
	   __eStateSolver_dirichlet_velocity,
	   __eStateSolver_dirichlet_marker,
	   __eStateSolver_dirichlet_pressure_symbolic,
	   __eStateSolver_dirichlet_velocity_symbolic,
	   __eStateSolver_dirichlet_marker_symbolic,
	   __eStateSolver_initial_condition,
	   __eStateSolver_slip_velocity_symbolic,
	   __eStateSolver_initial_condition_dg,
	   __eStateSolver_dirichlet_dg,
	   __eStateSolver_exit);
					  
#endif
