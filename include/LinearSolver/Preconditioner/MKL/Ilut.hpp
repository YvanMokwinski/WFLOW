#ifndef __header_LinearSolver_Preconditioner_MKL_Ilut_hpp__
#define __header_LinearSolver_Preconditioner_MKL_Ilut_hpp__

#include "LinearSolver/ILinearOperator.hpp"
#include "LinearSolver/Iterative/MKL/Parameters.hpp"

#define MKL_ILP64 1
#include "mkl.h"
#include "mkl_rci.h"

namespace LinearSolver
{
  namespace Preconditioner
  {
    namespace MKL
    {
  
      /// <summary>
      /// Implementation of the Incomplete LU Factorization with level 0 of fill-in.
      /// </summary>
      /// <remarks>
      /// Avoid using this preconditioner with the Conjugate Gradient sparse iterative solver because in general, 
      /// it produces a non-symmetric resulting matrix even if the original matrix is symmetric.
      /// </remarks>
      class Ilut : public IInverseOperator
      {

      private:
	struct Error
	{
	  /// <summary>
	  /// Enumerate error values.
	  /// </summary>
	public: typedef enum Value
	  {
	    /// <summary>
	    /// Indicates that the task completed normally.
	    /// </summary>
	    Completed = 0,
		   
	    /// <summary>
	    /// Indicates that the routine was interrupted because of an error: the number of elements 
	    /// in some matrix row specified in the sparse format is equal to or less than 0.
	    /// </summary>
	    EncounteredAnEmptyRow = -101,
		   
	    /// <summary>
	    /// Indicates that the routine was interrupted because the value of the computed diagonal 
	    /// element is less than the product of the given tolerance and the current matrix row norm, and it cannot be replaced as ipar[30]=0.
	    /// </summary>
	    InvalidDiagonalValue = -102,

	    /// <summary>
	    /// Indicates that the routine was interrupted because the element ia[i] is less than or equal to the element ia[i - 1] (see Sparse Matrix Storage Format).
	    /// </summary>
	    EncounteredAnInvalidOffset = -103,

	    /// <summary>
	    /// Indicates that the routine was interrupted because the memory is insufficient for the internal work arrays
	    /// </summary>
	    NotEnoughInternalMemory = -104,

	    /// <summary>
	    /// Indicates that the routine was interrupted because the input value of maxfil is less than 0.
	    /// </summary>
	    InvalidMaxfil = -105,

	    /// <summary>
	    /// Indicates that the routine was interrupted because the size n of the input matrix is less than 0.
	    /// </summary>
	    InvalidMatrixDimension = -106,

	    /// <summary>
	    /// Indicates that the routine was interrupted because an element of the array ja is less than 1, or greater than n (see Sparse Matrix Storage Format).
	    /// </summary>
	    InvalidMatrixIndexing = -107,

	    /// <summary>
	    /// The value of maxfil is greater than or equal to n. The calculation is performed with the value of maxfil set to (n-1).
	    /// </summary>
	    InvalidMaxfil2 = 101,

	    /// <summary>
	    /// The value of tol is less than 0. The calculation is performed with the value of the parameter set to (-tol)
	    /// </summary>
	    InvalidTolerance = 102,

	    /// <summary>
	    /// The absolute value of tol is greater than value of dpar[30]; it can result in instability of the calculation.
	    /// </summary>
	    IncompatibleToleranceAndDiagonalElementThreshold = 103,

	    /// <summary>
	    /// The value of dpar[30] is equal to 0. It can cause calculations to fail.
	    /// </summary>
	    InvalidDiagonalElementThreshold = 104,
	  } Type;

	  /// <summary>
	  /// Get the error message from the routine <see cref="ILUT"/>.
	  /// </summary>
	  /// <param name="error_">The parameter 'ierr' from the routine.</param>
	  /// <returns>The error message related to the error.</returns>
	public: static string GetMessage(const Error::Type error_)
	  {
	    switch (error_)
	      {
	      case Error::Completed:
		{
		  return
		    "Indicates that the task completed normally.";
		}

	      case Error::EncounteredAnEmptyRow:
		{
		  return
		    "Indicates that the routine was interrupted because of an error: the number of elements in some matrix row specified in the sparse format is equal to or less than 0.";
		}

	      case Error::InvalidDiagonalValue:
		{
		  return "Indicates that the routine was interrupted because the value of the computed diagonal element is less than the product of the given tolerance and the current matrix row norm, and it cannot be replaced as ipar[30]=0.";
		}

	      case Error::EncounteredAnInvalidOffset:
		{
		  return
		    "Indicates that the routine was interrupted because the element ia[i] is less than or equal to the element ia[i - 1] (see Sparse Matrix Storage Format)";
		}

	      case Error::NotEnoughInternalMemory:
		{
		  return
		    "Indicates that the routine was interrupted because the memory is insufficient for the internal work arrays.";
		}


	      case Error::InvalidMaxfil:
		{
		  return
		    "Indicates that the routine was interrupted because the input value of maxfil is less than 0.";
		}


	      case Error::InvalidMatrixDimension:
		{
		  return
		    "Indicates that the routine was interrupted because the size n of the input matrix is less than 0.";
		}


	      case Error::InvalidMatrixIndexing:
		{
		  return
		    "Indicates that the routine was interrupted because an element of the array ja is less than 1, or greater than n (see Sparse Matrix Storage Format).";
		}


	      case Error::InvalidMaxfil2:
		{
		  return
		    "The value of maxfil is greater than or equal to n. The calculation is performed with the value of maxfil set to (n-1).";
		}

	      case Error::InvalidTolerance:
		{
		  return
		    "The value of tol is less than 0. The calculation is performed with the value of the parameter set to (-tol)";
		}

	      case Error::IncompatibleToleranceAndDiagonalElementThreshold:
		{
		  return
		    "The absolute value of tol is greater than value of dpar[30]; it can result in instability of the calculation.";
		}


	      case Error::InvalidDiagonalElementThreshold:
		{
		  return "The value of dpar[30] is equal to 0. It can cause calculations to fail.";
		}


	      default:
		{
		  return "Unknown error value";
		}
	      }
	  };

	};
	
	/// <summary>
	/// The matrix A.
	/// </summary>
      private: pSparse m_A;
	
	/// <summary>
	/// The preconditioner.
	/// </summary>
      private: pSparse m_P;
	
	/// <summary>
	/// Tolerance for threshold criterion for the resulting entries of the preconditioner.
	/// </summary>
      private: R m_tol;
	
	/// <summary>
	/// Maximum fill-in, which is half of the preconditioner bandwidth. 
	/// The number of non-zero elements in the rows of the preconditioner cannot exceed (2*<see cref="m_maxfil"/>+1)
	/// </summary>
      private: I m_maxfil;
	
	/// <summary>
	/// Parameter from RCI FGMRES, first initialized by the FGMRES method.
	/// </summary>
      private: LinearSolver::Iterative::MKL::Parameters * m_parameters;
	  
	/// <summary>
	/// Status.
	/// </summary>
      private: Error::Type m_error;
	
	/// <summary>
	/// Get the required of number of coefficients of the incomplete LU factorization.
	/// </summary>
	/// <param name="n_">The dimension.</param>
	/// <param name="maxfil_">The maximum fill.</param>
	/// <returns>The required number of coefficients.</returns>
      private: static I GetRequiredNumberOfCoefficients(const I n_,
							const I maxfil_)
	{
	  return (2 * maxfil_ + 1) * n_ - maxfil_ * (maxfil_ + 1) + 1;
	};

	
	/// <see cref="IInverseOperator.SizeOfTemporaryVector"/>>
      public: I GetSizeOfTemporaryVector()const
	{
	  return Sparse_get_n(this->m_A);
	};

	
	/// <see cref="IInverseOperator.ErrorMessage"/>
      public: std::string GetErrorMessage()const
	{
	  return Error::GetMessage(m_error); 
	};
	
      public: virtual ~Ilut()
	{
	  this->m_A 		= NULL;
	  this->m_P 		= Sparse_kill(this->m_P);
	  this->m_parameters 	= NULL;
	  this->m_error 	= Error::Completed;	  
	};

	/// <summary>
	/// Create an Incomplete LU Factorization preconditioner.
	/// </summary>
	/// <param name="A">The matrix A from which we need to compute a preconditioner.</param>
	/// <param name="tol">Tolerance for threshold criterion for the resulting entries of the preconditioner.</param>
	/// <param name="maxfil">Maximum fill-in, which is half of the preconditioner bandwidth. 
	/// The number of non-zero elements in the rows of the preconditioner cannot exceed (2*<paramref name="maxfil"/>+1)</param>
	/// <param name="parameters">Parameters of the RCI FGMRES computations, first initialized by the MKL FGMRES.</param>
      public: Ilut(pSparse A_,
		   const R tol_,
		   const I maxfil_,
		   LinearSolver::Iterative::MKL::Parameters * parameters_)
	{
	  
	  //
	  // Assign data members
	  //
	  this->m_tol = tol_;
	  this->m_maxfil = maxfil_;
	  this->m_parameters = parameters_;
	  this->m_A = A_;

	  
	  pI paramIntegers 	= this->m_parameters->GetParamIntegers();
	  pR paramReals 	= this->m_parameters->GetParamReals();
	  
	  //
	  // ipar[30]
	  // specifies how the routine operates if the value of the computed diagonal element is less 
	  // than the current matrix row norm multiplied by the value of the parameter tol. 
	  // If ipar[30] = 0, then the calculation is stopped and the routine returns non-zero error value. 
	  // Otherwise, the value of the diagonal element is set to a value determined by dpar[30] 
	  // (see its description below), and the calculations continue.
	  //
	  paramIntegers[30] = 1;
	  
	  //
	  // used to adjust the value of small diagonal elements. 
	  // Diagonal elements with a value less than the current matrix row norm multiplied 
	  // by tol are replaced with the value of dpar[30] multiplied by the matrix row norm.
	  //
	  paramReals[30] = 1.0e-5;

	  //
	  // Get the dimension
	  //
	  const I N = Sparse_get_n(this->m_A);
	  
	  //
	  // Compute the required number of coefficient.
	  //
	  const I requiredNumberOfCoefficients = Ilut::GetRequiredNumberOfCoefficients(N,
										       this->m_maxfil);
	  

	  //
	  // The preconditioner P needs to be 1-based indexed.
	  //
	  this->m_P = Sparse_new(N,
				 N,
				 requiredNumberOfCoefficients);
	  this->m_P->format = 1;

	};

	/// <see cref="IInverseOperator.Compute"/>
      public: void Compute(bool* outHasFailed_)
	{
	  *outHasFailed_ = false;
	  this->m_error = Error::Completed;
	    
	  //
	  // Set the matrix to 1-based indexing.
	  //
	  const L isFortranIndexed = Sparse_is_fortran_index(m_A);
	  if (__emnsNO == isFortranIndexed)
	    {
	      Sparse_fortran_indexation(this->m_A);
	    }
	    
	  //
	  // Compute the incomplete factorization
	  //
	    
	  {
	    I n = Sparse_get_n(this->m_A);
	    I lierr;
	    dcsrilut(&n,
		     Sparse_get_ownx(this->m_A),
		     Sparse_get_ownb(this->m_A),
		     Sparse_get_owni(this->m_A),
		     Sparse_get_ownx(this->m_P),
		     Sparse_get_ownb(this->m_P),
		     Sparse_get_owni(this->m_P),
		     &this->m_tol,
		     &this->m_maxfil,
		     m_parameters->GetParamIntegers(),
		     m_parameters->GetParamReals(),
		     &lierr);
	    this->m_error = (Error::Type)lierr;
	  }
	    
	  //
	  // Re-set the matrix to 0-based indexing.
	  //
	  if (__emnsNO == isFortranIndexed)
	    {
	      Sparse_c_indexation(this->m_A);
	    }
	    
	  switch (this->m_error)
	    {
	    case Error::Completed:
	      {
		break;
	      }
	    default:
	      {
		*outHasFailed_ = true;
		break;
	      }
	    }
	    
	};

	/// <see cref="IInverseOperator.Apply"/>
      public: void Apply(const char*	transpose_,
			 pR 		y_,
			 cst_pR 	x_,
			 const I  	tmpSize_,
			 pR  		tmp_,
			 bool* 		outHasFailed_)
	{
      
	  *outHasFailed_ = false;
	  I dimension = Sparse_get_n(this->m_A);

	  const char * sL = "L";
	  const char * sU = "U";
	  const char * sN = "N";

	  //
	  // Forward substitution
	  //
	  mkl_dcsrtrsv(sL,
		       sN,
		       sU,
		       &dimension,
		       Sparse_get_ownx(this->m_P),
		       Sparse_get_ownb(this->m_P),
		       Sparse_get_owni(this->m_P),
		       x_,
		       tmp_);
	  
	  
	  //
	  // Backward substitution
	  //
	  mkl_dcsrtrsv(sU,
		       sN,
		       sN,
		       &dimension,
		       Sparse_get_ownx(this->m_P),
		       Sparse_get_ownb(this->m_P),
		       Sparse_get_owni(this->m_P),			   
		       tmp_,
		       y_);
	    
	};

	  
      };

    }; // namespace MKL

  }; // namespace Preconditioner
  
}; // namespace LinearSolver

#endif
