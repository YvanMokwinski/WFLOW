
#pragma once

struct SolverConstants
{
  
public:
  
  typedef R ValueType;
  
  typedef enum Enum
    {
      ERROR=0,
      iweber,
      ireynold,
      ifroude,
      vmin,
      vmax,
      coeff_ratio_viscosity,
      coeff_ratio_density,
      dt,
      dti,
      idt,
      midt,
      ratio_dti,
      iratio_dti,
      eps,
      ieps,      
      ALL_VALUES
    } EnumType;

  inline ValueType& operator[](const EnumType enumType)
  {
    return this->m_values[enumType];
  };

  inline const ValueType& operator[](const EnumType enumType) const
  {
    return this->m_values[enumType];
  };
  
protected:

  ValueType m_values[ALL_VALUES];
};
