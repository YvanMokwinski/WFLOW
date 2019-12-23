#include "eTensionMethod.h"

ConfigEnumStrings(eTensionMethod,
		  "no",
		  "csf",
		  "ssf");

ConfigOptionTable(L,
		  eTensionMethod,
		  {__eTensionMethod_NO,__emnsYES,"no tension model","--tension-no"},
		  {__eTensionMethod_CSF,__emnsNO,"CSF tension model","--tension-csf"},
		  {__eTensionMethod_SSF,__emnsNO,"SSF tension model","--tension-ssf"} );



