﻿#ifndef __LinearSolver_Preconditioner_Jacobi_hpp__
#define __LinearSolver_Preconditioner_Jacobi_hpp__

#include "LinearSolver/DiagonalMatrix.hpp"

namespace LinearSolver
{
  namespace Preconditioner
  {
    
    /// <summary>
    /// Implementation of the Jacobi preconditioner.
    /// </summary>
    /// <remarks>
    /// For a matrix A, the Jacobi preconditioner P is chosen to be the diagonal of matrix A: P:=Diagonal(A).
    /// The application of the Jacobi preconditioner is y:= P^{-1} x, which is straightforward since P is diagonal.
    /// </remarks>
    class Jacobi : public DiagonalMatrix, public IInverseOperator
    {

    private:
      /// <summary>
      /// The matrix A.
      /// </summary>
      pSparse m_A;

      /// <summary>
      /// The minimum of the absolute value when computing the inverse.
      /// </summary>
      R m_minAbsoluteValue;

    public:
      /// <summary>
      /// Constructor.
      /// </summary>
      /// <param name="A">The sparse matrix.</param>
      /// <param name="minAbsoluteValue">The minimum of the absolute value when computing the inverse.</param>
      Jacobi(pSparse A, 
	     const R minAbsoluteValue = ((R)1.1e-16))
	: DiagonalMatrix(NULL != A ? Sparse_get_n(A) : 1)
      {
	this->m_A = A;
	this->m_minAbsoluteValue = minAbsoluteValue;
      };
     

      /// <see cref="IInverseOperator.Compute"/>
      void Compute(bool* outHasFailed)
      {
	//
	// Basically, it cannot fail.
	//
	*outHasFailed = false;

	//
	// Extract the diagonal.
	//
	Sparse_extractDiagonal(this->m_A,
			       this->m_values);

	
	//
	// Compute the inverse.
	//
	Inverse(this->m_minAbsoluteValue);
      };
     
      /// <see cref="IInverseOperator.SizeOfTemporaryVector"/>>
      I GetSizeOfTemporaryVector()const
      {
	return 1;
      };
     
      /// <see cref="IInverseOperator.Apply"/>>
      void Apply(const char*	transpose,
		 pR 		y,
		 cst_pR 	x,
		 const I	tmpSize,
		 pR  		tmp,
		 bool* 		outHasFailed)
      {
       
	//
	// We do not care about the temporary vector.
	//
	*outHasFailed = false;
	this->Amux(transpose,
		   y,
		   x);   
      };
     
      /// <see cref="IInverseOperator.ErrorMessage"/>>
      std::string GetErrorMessage() const
      {
	return std::string("No error message available, this preconditioner should never fail.");
      };

    };

  };

};

#endif
