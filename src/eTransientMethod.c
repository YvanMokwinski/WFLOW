#include "eTransientMethod.h"

ConfigEnumStrings(eTransientMethod,
		  "no",
		  "euler",
		  "imr",
		  "gear_euler",
		  "gear_imr",
		  "trapeze");

ConfigOptionTable(L,
		  eTransientMethod,
		  {__eTransientMethod_NO,__emnsYES,"no tension model","--transient-no"},
		  {__eTransientMethod_EULER,__emnsNO,"no tension model","--transient-euler"},
		  {__eTransientMethod_IMR,__emnsNO,"no tension model","--transient-imr"},
		  {__eTransientMethod_GEAREULER,__emnsNO,"no tension model","--transient-gear-euler"},
		  {__eTransientMethod_GEARIMR,__emnsNO,"CSF tension model","--transient-gear-imr"},
		  {__eTransientMethod_TRAPEZE,__emnsNO,"SSF tension model","--transient-trapeze"} );



