#include "eFileSuffix.h"

ConfigEnumStrings(eFileSuffix,
		  ".cmgdata",
		  ".mnsbin",
		  ".tex",
		  ".m",
		  ".mesh",
		  ".bb",
		  ".pie",
		  ".txt",
		  ".vtk",
		  ".mnsbin.gz",
		  ".mesh.gz",
		  ".bb.gz");




cst_eFileSuffix eFileSuffix_get_suffix			(cst_pS const 	filename_)
{
  eFileSuffix fileSuffix = __eFileSuffix_ERROR;
  I lenFilename = 0;
  I lenSuffix = 0;
  while (filename_[lenFilename]!='\0') ++lenFilename;
  
  for (++fileSuffix;fileSuffix<__eFileSuffix_ALL;++fileSuffix)
    {
      cst_pS suffix = eFileSuffix_get_string(fileSuffix);

      lenSuffix = 0;
      while (suffix[lenSuffix]!='\0') ++lenSuffix;
      if (!strcmp(suffix,&filename_[lenFilename-lenSuffix]))
	{
	  break;
	}	
    } 
  if (fileSuffix<__eFileSuffix_ALL)
    {	
      return fileSuffix;
    }
  else
    {
      return __eFileSuffix_ERROR;
    } 

}

cst_eFileSuffix eFileSuffix_get_suffixAndBasename		(cst_pS const 	filename_,
								 STR 		basename_)
{
  eFileSuffix fileSuffix = __eFileSuffix_ERROR;
  I lenFilename = 0;
  I lenSuffix = 0;
  while (filename_[lenFilename]!='\0') ++lenFilename;
  
  for (++fileSuffix;fileSuffix<__eFileSuffix_ALL;++fileSuffix)
    {
      cst_pS suffix = eFileSuffix_get_string(fileSuffix);

      lenSuffix = 0;
      while (suffix[lenSuffix]!='\0') ++lenSuffix;
      if (!strcmp(suffix,&filename_[lenFilename-lenSuffix]))
	{
	  break;
	}	
    } 
  if (fileSuffix<__eFileSuffix_ALL)
    {	
      sprintf(basename_,"%s",filename_);

      
      switch(fileSuffix)
	{
	case   __eFileSuffix_CmgData:
	case   __eFileSuffix_MnsBinary:
	case   __eFileSuffix_Latex:
	case   __eFileSuffix_Matlab:
	case   __eFileSuffix_MeditMesh:
	case   __eFileSuffix_MeditSolution:
	case   __eFileSuffix_Pirate:
	case   __eFileSuffix_Text:
	case   __eFileSuffix_VtkAscii:
	  {
	    basename_[lenFilename - lenSuffix] = '\0';
	    break;
	  }
	case   __eFileSuffix_ZippedMnsBinary:
	case   __eFileSuffix_ZippedMeditMesh:
	case   __eFileSuffix_ZippedMeditSolution:
	  {
	    basename_[lenFilename - 3] = '\0';
	    break;
	  }
	case   __eFileSuffix_ERROR:
	case   __eFileSuffix_ALL:
	  {
	    break;
	  }
	}

      return fileSuffix;
    }
  else
    {
      sprintf(basename_,"%s",filename_);
      return __eFileSuffix_ERROR;
    } 
}
