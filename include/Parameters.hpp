#ifndef __header_Parameters_hpp__
#define __header_Parameters_hpp__

#include "Cmdline.h"

#include "ns_enum.h"
#include "eHeaviside.h"
#include "eDirac.h"


template <typename T,typename V> struct OptionBase
{
public:
  T flag;
  V default_value;
  const char * description;
  const char * name_cmdline;
  const char * name_configfile;
};




struct KindTransientMethod
{
public:
  typedef enum Enum
    {
      ERROR=0,
      UNDEFINED,
      EULER,
      IMR,
      GEAREULER,
      GEARIMR,
      TRAPEZE,
      ALL_VALUES
    } EnumType;
};

  
struct KindTensionMethod
{
public:
  typedef enum Enum
    {
      ERROR=0,
      UNDEFINED,
      CSF,
      SSF,
      ALL_VALUES
    } EnumType;
};


struct KindEquation
{
public:
  typedef enum Enum
    {
      ERROR=0,
      VELOCITY,
      PRESSURE,
      TRANSPORT,
      ALL_VALUES
    } EnumType;
};

struct KindLinearSolver
{
public:
  typedef enum Enum
    {
      ERROR=0,
      Direct,
      Iterative,
      ALL_VALUES
    } EnumType;
};

struct InfoString
{
public:
  typedef const char * ValueType;
  typedef enum Enum
    {
      ERROR=0,     
      name,
      pblmname,
      ofilename,
      ALL_VALUES
    } EnumType;

  inline static const OptionBase<EnumType,ValueType>*GetTable()
  {
    static const OptionBase<EnumType,ValueType> s_opts[ALL_VALUES]=
      {
	{ERROR,NULL,"","",""},
	{name,"widji","name of the process","--name","NAME"},
	{pblmname,"laplace","name of the problem","--pblm","PBLMNAME"},
	{ofilename,"ns.out","name of the output file","-o","OFILENAME"}
      };
    return s_opts;
  };
  
};


/// <summary>
/// Enumerate the integer info.
/// </summary>
struct InfoInteger
{
public:
  typedef I ValueType;
  typedef enum Enum
    {
      ERROR=0,
      newton_maxiter,
      ntime,
      noboundary_vcod,
      mesh_adaptivity_maxiter,
      ntime_interval,
      heaviside,
      dirac,
      restart,
      nproc,
      ALL_VALUES
    } EnumType;

  inline static const OptionBase<EnumType,ValueType>*GetTable()
  {
    static const OptionBase<EnumType,ValueType> s_opts[ALL_VALUES]=
      {
	{ERROR,0,"","",""},
	{newton_maxiter,((I)20),"newton max iter","--newton-niter","NEWTON_MAXITER"},
	{ntime,((I)7),"n time step","--ntime","NTIME"},
	{noboundary_vcod,((I)100),"no boundary cod","--noboundary-cod","NOBOUNDARY_COD"},
	{mesh_adaptivity_maxiter,((I)5),"n time step","--adapt-niter","ADAPT_NITER"},
	{ntime_interval,((I)5),"ntime interval","--ns-ntime-interval","NTIME_INTERVAL"},
	{heaviside,((I)__eHeaviside_m0p6),"heaviside smooth","--ns-heaviside","HEAVISIDE"},
	{dirac,((I)__eDirac_m0p6),"dirac smooth","--ns-dirac","DIRAC"},
	{restart,((I)0),"restart time step","--ns-restart","RESTART"},
	{nproc,((I)1),"nproc","-n","NPROC"}
      };
    return s_opts;
  };

  
  
};

/// <summary>
/// Enumerate the logical info.
/// </summary>
struct InfoReal
{
public:

  typedef double ValueType;
  typedef enum Enum
    {
      ERROR=0,
      hmin,
      hmax,
      ratio_viscosity,
      ratio_density,
      dtmin,
      dtmax,
      epsmin,
      epsmax,
      newton_tol_residu,
      newton_tol_correc,
      reynold,
      weber,
      froude,
      ALL_VALUES
    } EnumType;

  inline static const OptionBase<EnumType,ValueType>*GetTable()
  {
    static const OptionBase<EnumType,ValueType> s_opts[ALL_VALUES]=
      {
	{ERROR,((R)0.0),"","",""},
	{hmin,((R)1.0e-4),"hmin","--hmin","HMIN"},
	{hmax,((R)0.25),"hmax","--hmax","HMAX"},
	{ratio_viscosity,((R)1.0),"ratio viscosity","--nu","NU"},
	{ratio_density,((R)1.0),"ratio density","--rho","RHO"},
	{dtmin,((R)1.0e-4),"dtmin","--dt","DT"},
	{dtmax,((R)0.5),"dtmax","--dtmax","DTMAX"},
	{epsmin,((R)0.25),"epsmin","--epsmin","EPSMIN"},
	{epsmax,((R)0.9),"epsmax","--epsmax","EPSMAX"},
	{newton_tol_residu,((R)1.0e-8),"stopping criteria based on residu","--ntolr","NEWTON_TOLR"},
	{newton_tol_correc,((R)1.0e-5),"stopping criteria based on correction","--ntolc","NEWTON_TOLC"},
	{reynold,((R)0.5),"reynold number","--reynold","REYNOLD"},
	{weber,((R)0.5),"weber number","--weber","WEBER"},
	{froude,((R)0.5),"froude number","--froude","FROUDE"}
      };
    return s_opts;
  };
};

/// <summary>
/// Enumerate the logical info.
/// </summary>
struct InfoLogical
{
  
public:
  
  typedef bool ValueType;
  typedef enum Enum
    {
      ERROR=0,
      verbose,
      mesh_adaptivity,
      time_adaptivity,
      transport,
      transport_uncoupled,
      transport_galerkin,
      pressure_uncoupled,
      pressure_freematrix,
      slip,
      color,
      tension,
      redistance,
      pspg,
      usupg,
      axisymetric_x,
      axisymetric_y,
      dynamic_boundary_condition,
      vnsdg,
      skip_verif,
      ALL_VALUES
    } EnumType;
  
  inline static const OptionBase<EnumType,bool>*GetTable()
  {
    static const OptionBase<EnumType,bool> s_opts[ALL_VALUES]=
      {
	{ERROR,false,"","",""},
	{verbose,false,"activate verbose","-v","VERBOSE"},
	{mesh_adaptivity,false,"mesh adaptivity","--mesh-adaptivity","MESH_ADAPTIVITY"},
	{time_adaptivity,false,"time adaptivity","--time-adaptivity","TIME_ADAPTIVITY"},
	{transport,true,"apply transport","--transport","TRANSPORT"},
	{transport_uncoupled,false,"apply uncoupled transport","--transport-uncoupled","TRANSPORT_UNCOUPLED"},
	{transport_galerkin,true,"apply dg for transport","--poulou","DG"},
	{pressure_uncoupled,false,"uncouple pressure","--pressure-uncoupled","PRESSURE_UNCOUPLED"},
	{pressure_freematrix,false,"do not build pressure B matrix ","--pfreem","PFREEM"},
	{slip,false,"apply slip boundary condition","--slip","SLIP"},
	{color,true,"choose color signed distance","--sd","SD"},
	{tension,true,"do not apply tension","--notension","NOTENSION"},
	{redistance,false,"apply redistance algorithm","--redistance","REDISTANCE"},
	{pspg,false,"apply PSPG to pressure equations","--pspg","PSPG"},
	{usupg,false,"apply SUPG to velocity equations","--usupg","USUPG"},
	{axisymetric_x,false,"set axisymetric x","--axix","AXIX"},
	{axisymetric_y,false,"set axisymetric y","--axiy","AXIY"},
	{dynamic_boundary_condition,false,"no no ","--no","DYNAMIC_BOUNDARY_CONDITION"},
	{vnsdg,false,"enable vnsdg","--dg","VNSDG"},
	{skip_verif,false,"skip verif","-f","SKIP_VERIF"}
      };
    return s_opts;
  };
  
};


inline InfoLogical::EnumType & operator ++(InfoLogical::EnumType & self_) { self_=static_cast<InfoLogical::EnumType>(self_+1); return self_; };
inline InfoReal::EnumType & operator ++(InfoReal::EnumType & self_) { self_=static_cast<InfoReal::EnumType>(self_+1); return self_; };
inline InfoInteger::EnumType & operator ++(InfoInteger::EnumType & self_) { self_=static_cast<InfoInteger::EnumType>(self_+1); return self_; };
inline InfoString::EnumType & operator ++(InfoString::EnumType & self_) { self_=static_cast<InfoString::EnumType>(self_+1); return self_; };



template <typename _Info> struct InfoAccess
{
public:
  typedef typename _Info::EnumType EnumType;
  typedef typename _Info::ValueType ValueType;
  
  static const char * GetNameCmdline(const EnumType value_)
  {
    static const OptionBase<EnumType,ValueType>* s_table = _Info::GetTable();
    return s_table[value_].name_cmdline;
  };

  static ValueType GetDefaultValue(const EnumType value_)
  {
    static const OptionBase<EnumType,ValueType>* s_table = _Info::GetTable();
    return s_table[value_].default_value;
  };


};


typedef InfoAccess<InfoLogical> InfoAccessLogical;
typedef InfoAccess<InfoReal> InfoAccessReal;
typedef InfoAccess<InfoInteger> InfoAccessInteger;
typedef InfoAccess<InfoString> InfoAccessString;

class Parameters
{
private:
  bool 				m_linfo[InfoLogical::ALL_VALUES];
  I 				m_iinfo[InfoInteger::ALL_VALUES];
  R				m_rinfo[InfoReal::ALL_VALUES];
  std::string			m_sinfo[InfoString::ALL_VALUES];
  KindTransientMethod::EnumType	m_method_transient_scheme[KindEquation::ALL_VALUES];
  KindLinearSolver::EnumType	m_method_linearsolver[KindEquation::ALL_VALUES];
  KindTensionMethod::EnumType	m_tension_method;

public:

  Parameters(pCmdline const cmdline_,const bool	have_configfile_)
  {

    { InfoString::EnumType i = InfoString::ERROR;
      for (++i;i<InfoString::ALL_VALUES;++i)
	{
#if 0
	  if (__tableoptions_ens_sinfo[i].flag!=i)
	    {
	      fprintf(stderr,"ens_sinfo_argv:wrong dat nbase\n");
	      
	    }
#endif
	  STR name;
	  const L a = Cmdline_get_string(cmdline_,
					 InfoAccessString::GetNameCmdline(i),
					 &name[0]);
	  if ( (NOT a) AND (NOT have_configfile_) )
	    {
	      this->m_sinfo[i] = InfoAccessString::GetDefaultValue(i);
	    }
	} }


    { InfoInteger::EnumType i = InfoInteger::ERROR;
    for (++i;i<InfoInteger::ALL_VALUES;++i)
      {
#if 0
	if (__tableoptions_ens_rinfo[i].flag!=i)
	  {
	    fprintf(stderr,"ens_rinfo_argv:wrong dat nbase\n");
	  }
#endif
	const L a = Cmdline_get_integer(cmdline_,
					InfoAccessInteger::GetNameCmdline(i),
					&this->m_iinfo[i]);
	if ( (NOT a) AND (NOT have_configfile_) )
	  {
	    this->m_iinfo[i] = InfoAccessInteger::GetDefaultValue(i);
	  }
      } }

    
    { InfoReal::EnumType i = InfoReal::ERROR;
      for (++i;i<InfoReal::ALL_VALUES;++i)
	{
#if 0
	  if (__tableoptions_ens_rinfo[i].flag!=i)
	    {
	      fprintf(stderr,"ens_rinfo_argv:wrong dat nbase\n");
	    }
#endif
	  const L a = Cmdline_get_real(cmdline_,
				       InfoAccessReal::GetNameCmdline(i),
				       &this->m_rinfo[i]);
	  if ( (NOT a) AND (NOT have_configfile_) )
	    {
	      this->m_rinfo[i] = InfoAccessReal::GetDefaultValue(i);
	    }
	} }
    
    
    { InfoLogical::EnumType value = InfoLogical::ERROR;
      for (++value;value<InfoLogical::ALL_VALUES;++value)
	{
#if 0
	  if (__tableoptions_ens_linfo[i].flag!=i)
	    {
	      fprintf(stderr,"ens_linfo_argv:wrong dat nbase\n");
	    }
#endif	
	  const L a = Cmdline_get_logical(cmdline_,
					  InfoAccessLogical::GetNameCmdline(value));
	  if (a)
	    {
	      this->m_linfo[value] = !InfoAccessLogical::GetDefaultValue(value);
	    }
	  else
	    {
	      if ( (NOT a) AND (NOT have_configfile_) )
		{
		  this->m_linfo[value] = InfoAccessLogical::GetDefaultValue(value);
		}
	    }
	} }

#if 0
    ens_linfo_from_Cmdline(this->m_linfo,have_configfile_ ? __emnsNO : __emnsYES,cmdline_);
    ens_iinfo_from_Cmdline(this->m_iinfo,have_configfile_ ? __emnsNO : __emnsYES,cmdline_);
    ens_rinfo_from_Cmdline(this->m_rinfo,have_configfile_ ? __emnsNO : __emnsYES,cmdline_);
    ens_sinfo_from_Cmdline(this->m_sinfo,have_configfile_ ? __emnsNO : __emnsYES,cmdline_);
#endif
  };

  virtual ~Parameters()
  {
  };

  inline bool GetInfoLogical(const InfoLogical::EnumType value_) const
  {
    return this->m_linfo[value_];
  };

  inline I GetInfoInteger(const InfoInteger::EnumType value_) const
  {
    return this->m_iinfo[value_];
  };
  
  inline R GetInfoReal(const InfoReal::EnumType value_) const
  {
    return this->m_rinfo[value_];
  };
  
  inline void SetInfoReal(const InfoReal::EnumType infoReal_,const R value_)
  {
    this->m_rinfo[infoReal_] = value_;
  };

  inline void SetInfoInteger(const InfoInteger::EnumType info_,const I value_)
  {
    this->m_iinfo[info_] = value_;
  };

  inline void SetInfoLogical(const InfoLogical::EnumType info_,const bool value_)
  {
    this->m_linfo[info_] = value_;
  };
  
  inline std::string GetInfoString(const InfoString::EnumType value_) const
  {
    return this->m_sinfo[value_];
  };

  inline void SetInfoString(const InfoString::EnumType info_,const std::string& value_)
  {
    this->m_sinfo[info_] = value_;
  };
  
  

  
  inline void SetKindLinearSolver(const KindEquation::EnumType kindEquation_,KindLinearSolver::EnumType value_) 
  {
    this->m_method_linearsolver[kindEquation_] = value_;
  };

  inline KindLinearSolver::EnumType GetKindLinearSolver(const KindEquation::EnumType kindEquation_) const
  {
    return this->m_method_linearsolver[kindEquation_];
  };

  inline KindTransientMethod::EnumType GetKindTransientMethod(const KindEquation::EnumType kindEquation_) const
  {
    return this->m_method_transient_scheme[kindEquation_];
  };

  inline KindTensionMethod::EnumType GetKindTensionMethod() const
  {
    return this->m_tension_method;
  };

  
  inline  void SetKindTransientMethod(const KindEquation::EnumType kindEquation_,const KindTransientMethod::EnumType value_) 
  {
     this->m_method_transient_scheme[kindEquation_] = value_;
  };

  inline  void SetKindTensionMethod(const KindTensionMethod::EnumType value_) 
  {
     this->m_tension_method = value_;
  };


  
};

extern "C"
{

  typedef void * pParameters;
  typedef const void * pParametersReadOnly;
  
  L 			Parameters_getl			(const void* 	self_,
							 const ens_linfo linfo_);
  R 			Parameters_getr			(const void* 	self_,
							 const ens_rinfo rinfo_);
  I 			Parameters_geti			(const void* 	self_,
							 const ens_iinfo iinfo_);
  
  eLinearSolver 	Parameters_get_eLinearSolver	(const void*self_,
							 eKindEquation	kindEquation_);
  eTransientMethod 	Parameters_get_eTransientScheme	(const void*self_,
							 eKindEquation	kindEquation_);

}


#endif
