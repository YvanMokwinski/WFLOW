#ifndef __header_IInverseOperator_hpp__
#define __header_IInverseOperator_hpp__

namespace LinearSolver
{

  /// \brief Abstract base class for all solvers.
  /// An InverseOperator computes the solution of \f$ A(x)=b\f$ where                                    
  /// \f$ A : X \to Y \f$ is an operator.
  /// Note that the solver "knows" which operator
  /// to invert and which preconditioner to apply (if any). The
  /// user is only interested in inverting the operator.
  /// InverseOperator might be a Newton scheme, a Krylov subspace method, or a direct solver or just anything.
  class IInverseOperator
  {
  public:
  
    /// <summary>
    /// The error message if any operation failed.
    /// </summary>
    virtual std::string GetErrorMessage() const = 0;
  
    /// <summary>
    /// The size of the required temporary array of double to apply the linear operator.
    /// </summary>
    /// <returns>The size of the required temporary vector to apply the linear operator.</returns>
    virtual I 	GetSizeOfTemporaryVector() const = 0;

  
    /// <summary>
    /// Compute the inverse operator.
    /// </summary>
    /// <param name="outHasFailed">Indicates if computing the inverse operator has failed. 
    /// True if it has failed, false otherwise.</param>
    virtual void 	Compute(bool* outHasFailed) = 0;
  
    /// <summary>
    /// Apply the inverse operator: y := Op(x).
    /// </summary>
    /// <param name="transpose">No transpose ('N' or 'n'), Transpose ('T' or 't').</param>
    /// <param name="y">The output vector.</param>
    /// <param name="x">The input vector.</param>
    /// <param name="tmpSize">The size of the temporary vector.</param>
    /// <param name="tmp">The temporary array of double.</param>
    /// <param name="outHasFailed">Indicates if applying the operator has failed. 
    /// True if it has failed, false otherwise.</param>
    virtual void 	Apply(const char*	transpose,
			      pR y,
			      cst_pR x,
			      const I  		tmpSize,
			      pR  		tmp,
			      bool* 		outHasFailed) = 0;
  
  };

};
  
#endif
