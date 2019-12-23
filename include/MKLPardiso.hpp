#ifndef __header_MKLPardiso_hpp__
#define __header_MKLPardiso_hpp__

#include "IInverseOperator.hpp"

/// <summary>
/// Implementation of the linear solver with PARDISO.
/// </summary>
class MKLPardiso : public IInverseOperator
{
private:

  /// <summary>
  /// Type of the matrix.
  /// </summary>
  static const Int32 s_mtype = (int) MKLWrapper.DirectSolvers.PARDISO.MatrixType.ReadAndUnsymmetric;
  
  /// <summary>
  /// Maximal number of factors in memory. Generally used value is 1.
  /// </summary>
  static const Int32 s_maxfct = 1;
  
  /// <summary>
  /// The number of matrix (from 1 to <see cref="s_maxfct"/>) to solve.
  /// </summary>
  static const Int32 s_mnum = 1;

  /// <summary>
  /// Message level information.
  /// </summary>
  /// <remarks>
  /// 0: PARDISO generates no output
  /// 1: PARDISO prints statistical information
  /// </remarks>
  static const Int32 s_msglvl = 0;

  /// <summary>
  /// We fix the number of right-hand sides.
  /// </summary>
  static const Int32 s_nrhs = 1;

  
protected:
  /// <summary>
  /// The error variable.
  /// </summary>
  MKLPARDISOError m_error;

  /// <summary>
  /// Array of pointers.
  /// </summary>
  MKLPARDISO.Parameters m_parameters = new MKLPARDISO.Parameters();
  
  /// <summary>
  /// The matrix.
  /// </summary>
  SparseMatrixCSR32 m_A;

public:
  
  /// <summary>
  /// Constructor.
  /// </summary>
  /// <param name="A">The sparse matrix  <see cref="SparseMatrixCSR32"/>.</param>
  MKLPardiso(SparseMatrixCSR32 A)
  {


    m_A = A;
    m_error = MKLPARDISOError.NoError;

    //
    // Ensure we have C-indexing
    //
    m_parameters.FortranIndexing = false;

    //
    // Compute the symbolic factorization.
    // 
    DebugVerbose("symbolic factorization ...");

    MKLPARDISO.Pardiso(m_parameters,
		       maxfct: s_maxfct,
		       mnum: s_mnum,
		       mtype: s_mtype,
		       job: MKLPARDISOJob.SymbolicFactorization,
		       n: m_A.Dimension,
		       a: m_A.ValuesPtr,
		       ia: m_A.OffsetsPtr,
		       ja: m_A.IndicesPtr,
		       perm: IntPtr.Zero,
		       nrhs: s_nrhs,
		       msglvl: s_msglvl,
		       b: IntPtr.Zero,
		       x: IntPtr.Zero,
		       outError: out m_error);
    
    if (m_error != MKLPARDISOError.NoError)
      {
	throw new Exception(string.Format("MKL PARDISO has failed during symbolic factorization:{0}",
					  MKLPARDISO.GetErrorMessage(m_error)));
      }
    
    DebugVerbose("symbolic factorization done.");
  };

  
  /// <see cref="IInverseOperator.SizeOfTemporaryVector"/>>
  I GetSizeOfTemporaryVector() const
  {
    get { return 0; }
  };

  /// <see cref="IInverseOperator.Apply"/>
  virtual void Apply(char transpose,
		     void* y,
		     void* x,
		     Int32 tmpSize,
		     void* tmp,
		     bool* outHasFailed)
  {
    //
    // We do not care about the temporary vector tmp.
    //

    DebugVerbose("solve ...");

    m_parameters.Transpose = transpose == 'T' || transpose == 't';
    MKLPARDISO.Pardiso(m_parameters,
		       maxfct : s_maxfct,
		       mnum : s_mnum,
		       mtype : s_mtype,
		       job : MKLPARDISOJob.Solve,
		       n : m_A.Dimension,
		       a : m_A.ValuesPtr,
		       ia : m_A.OffsetsPtr,
		       ja : m_A.IndicesPtr,
		       perm : IntPtr.Zero,
		       nrhs : s_nrhs,
		       msglvl : s_msglvl,
		       b : x,
		       x : y,
		       outError : out m_error);
    outHasFailed = m_error != MKLPARDISOError.NoError;

    DebugVerbose("solve done.");
  };
  
  /// <see cref="IInverseOperator.Compute"/>
  void Compute(out bool outHasFailed)
  {
    MKLPARDISO.Pardiso(m_parameters,
		       maxfct : s_maxfct,
		       mnum : s_mnum,
		       mtype : s_mtype,
		       job : MKLPARDISOJob.NumericalFactorization,
		       n : m_A.Dimension,
		       a : m_A.ValuesPtr,
		       ia : m_A.OffsetsPtr,
		       ja : m_A.IndicesPtr,
		       perm : IntPtr.Zero,
		       nrhs : s_nrhs,
		       msglvl : s_msglvl,
		       b : IntPtr.Zero,
		       x : IntPtr.Zero,
		       outError : out m_error);

    outHasFailed = m_error != MKLPARDISOError.NoError;
    DebugVerbose("numerical factorization done.");
  }
  

  void GetErrorMessage(char message_[512])
  {
    MKLPARDISO.GetErrorMessage(m_error);
  }
  
  virtual ~MKLPardiso()
  {
    MKLPARDISO.Pardiso(m_parameters,
		       maxfct : s_maxfct,
		       mnum : s_mnum,
		       mtype : s_mtype,
		       job : MKLPARDISOJob.ReleaseMemory,
		       n : m_A.Dimension,
		       a : IntPtr.Zero,
		       ia : m_A.OffsetsPtr,
		       ja : m_A.IndicesPtr,
		       perm : IntPtr.Zero,
		       nrhs : s_nrhs,
		       msglvl : s_msglvl,
		       b : IntPtr.Zero,
		       x : IntPtr.Zero,
		       outError : out m_error);

  };


private:
  
  static void DebugVerbose(string message)
  {
    std::cout << "MKLPardiso: " << message << std::endl;
  };
  
  
};

#endif
