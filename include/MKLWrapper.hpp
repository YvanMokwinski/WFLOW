#ifndef __header_MKLWrapper_hpp__
#define __header_MKLWrapper_hpp__
   

/// <summary>
/// Available Preconditioners.
/// </summary>
class Preconditioners
{
  
  /// <summary>
  /// ILU0 preconditioner based on incomplete LU factorization of a sparse matrix.
  /// </summary>
  static class ILU0
  {
    /// <summary>
    /// Enumerate error values.
    /// </summary>
    enum Error
      {
	/// <summary>
	/// Indicates that the task completed normally.
	/// </summary>
	Completed = 0,
	
	/// <summary>
	/// Indicates that the routine was interrupted and that error occurred: 
	/// at least one diagonal element is omitted from the matrix in CSR3 format 
	/// </summary>
	MissingDiagonalElement = -101,
	
	/// <summary>
	/// Indicates that the routine was interrupted because the matrix contains a diagonal element with the value of zero.
	/// </summary>
	ZeroDiagonalElement = -102,
	
	/// <summary>
	/// Indicates that the routine was interrupted because the matrix contains a diagonal element which is so small that it could cause an overflow, 
	/// or that it would cause a bad approximation to ILU0.
	/// </summary>
	TooSmallDiagonalElement = -103,
	
	/// <summary>
	/// Indicates that the routine was interrupted because the memory is insufficient for the internal work array.
	/// </summary>
	NotEnoughInternalMemory = -104,
	
	/// <summary>
	/// Indicates that the routine was interrupted because the size n of the input matrix is less than 0.
               /// </summary>
	InvalidMatrixDimension = -105,
	
	/// <summary>
	/// Indicates that the routine was interrupted because the column indices ja are not in the ascending order.
	/// </summary>
	InvalidMatrixIndexing = -106,
      };
    
    /// <summary>
    /// Get the error message from the routine <see cref="ILUT"/>.
    /// </summary>
    /// <param name="error">The parameter 'ierr' from the routine.</param>
    /// <returns>The error message related to the error.</returns>
    static const char* ErrorMessage(Error error)
    {
      switch (error)
	{
	case Error.Completed:
	  {
	    return
	      "Indicates that the task completed normally.";
	  }
	  
	case Error.MissingDiagonalElement:
	  {
	    return
	      "Indicates that the routine was interrupted and that error occurred: at least one diagonal element is omitted from the matrix in CSR3 format.";
	  }
	  
	case Error.ZeroDiagonalElement:
	  {
	    return "Indicates that the routine was interrupted because the matrix contains a diagonal element with the value of zero.";
	  }
	  
	case Error.TooSmallDiagonalElement:
	  {
	    return
	      "Indicates that the routine was interrupted because the matrix contains a diagonal element which is so small that it could cause an overflow, or that it would cause a bad approximation to ILU0.";
	  }
	  
	case Error.NotEnoughInternalMemory:
	  {
	    return
	      "Indicates that the routine was interrupted because the memory is insufficient for the internal work arrays.";
	  }
	  
	  
	case Error.InvalidMatrixDimension:
	  {
	    return
	      "Indicates that the routine was interrupted because the size n of the input matrix is less than 0.";
	  }
	  
	  
	case Error.InvalidMatrixIndexing:
	  {
	    return
	      "Indicates that the routine was interrupted because the column indices ja are not in the ascending order.";
	  }
	  
	default:
	  {
	    return "Unknown error value";
	  }
	}
    };
    
    /// <summary>
    /// ILU0 preconditioner based on incomplete LU factorization of a sparse matrix.
    /// </summary>
    /// <remarks>
    /// The routine dcsrilu0 computes a preconditioner B of a given sparse matrix A stored in the format accepted in the direct sparse solvers:
    /// A~B=L*U , where L is a lower triangular matrix with a unit diagonal, U is an upper triangular matrix with a non-unit diagonal, and the portrait 
    /// of the original matrix A is used to store the incomplete factors L and U.
    /// CAUTION: This routine supports only one-based indexing of the array parameters.
    /// </remarks>
    /// 
    /// <param name="n">
    /// Size (number of rows or columns) of the original square  <paramref name="n"/>-by-<paramref name="n"/> matrix A.
    /// </param>
    /// 
    /// <param name="a">
    /// Array containing the set of elements of the matrix A. Its length is equal to the number of non-zero elements in the matrix A. 
    /// Refer to the values array description in the Sparse Matrix Storage Format for more details.
    /// </param>
    /// 
    /// <param name="ia">Array of size (<paramref name="n"/>+1) containing begin indices of rows of the matrix A such that <paramref name="ia"/>(i) is the index 
    /// in the array a of the first non-zero element from the row i. The value of the last element <paramref name="ia"/>(<paramref name="n"/>+1) is equal to the 
    /// number of non-zero elements in the matrix A, plus one. Refer to the rowIndex array description in the Sparse Matrix 
    /// Storage Format for more details.
    /// </param>
    /// 
    /// <param name="ja">Array containing the column indices for each non-zero element of the matrix A. 
    /// It is important that the indices are in increasing order per row. 
    /// The matrix size is equal to the size of the array <paramref name="a"/>. 
    /// Refer to the columns array description in the Sparse Matrix Storage Format for more details.
    /// CAUTION: If column indices are not stored in ascending order for each row of matrix, the result of the routine might not be correct.
    /// </param>
    /// 
    /// <param name="bilu0">
    /// Array B containing non-zero elements of the resulting preconditioning matrix B, stored in the format accepted in direct sparse solvers. 
    /// Its size is equal to the number of non-zero elements in the matrix A. Refer to the values array description in the Sparse Matrix Storage Format section for more details.
    /// </param>
    /// 
    /// <param name="ipar">
    /// Array of size 128. 
    /// This parameter specifies the integer set of data for both the ILU0 and RCI FGMRES computations. 
    /// Refer to the <paramref name="ipar"/> array description in the FGMRES Common Parameters for more details on FGMRES parameter entries. 
    /// The entries that are specific to ILU0 are listed below.
    /// 
    /// <paramref name="ipar"/> [30] specifies how the routine operates when a zero diagonal element occurs during calculation. 
    /// If this parameter is set to 0 (the default value set by the routine dfgmres_init), then the calculations are stopped and the routine returns a non-zero error value. 
    /// Otherwise, the diagonal element is set to the value of <paramref name="dpar"/>[31] and the calculations continue.
    /// 
    /// Note: You can declare the <paramref name="ipar"/>  array with a size of 32. However, for future compatibility you must declare the array <paramref name="ipar"/> with length 128
    /// </param>
    /// 
    /// <param name="dpar">Array of size 128. 
    /// This parameter specifies the double precision set of data for both the ILU0 and RCI FGMRES computations. 
    /// Refer to the <paramref name="dpar"/> array description in the FGMRES Common Parameters for more details on FGMRES parameter entries. 
    /// The entries specific to ILU0 are listed below:
    /// 
    /// <paramref name="dpar"/>[30] specifies a small value, which is compared with the computed diagonal elements. When <paramref name="ipar"/>[30] is not 0, 
    /// then diagonal elements less than <paramref name="dpar"/>[30] are set to <paramref name="dpar"/>[31]. The default value is 1.0e-16.
    /// Note: This parameter can be set to the negative value, because the calculation uses its absolute value.
    /// If this parameter is set to 0, the comparison with the diagonal element is not performed.
    /// 
    /// <paramref name="dpar"/>[31] specifies the value that is assigned to the diagonal element if its value is less than <paramref name="dpar"/>[30] (see above). The default value is 1.0e-10.
    /// </param>
    /// 
    /// <param name="outError">Error flag, gives information about the routine completion.</param>
     public static void Compute(Int32 n,
				IntPtr a,
				IntPtr ia,
				IntPtr ja,
				IntPtr bilu0,
				Int32[] ipar,
				double[] dpar,
				Error* outError)
    {
      int lierr;
      dcsrilu0(&n,
	       a,
	       ia,
	       ja,
	       bilu0,
	       ipar,
	       dpar,
	       &lierr);
      outError = (Error)lierr;               
    };
    
    
  };
  
  /// <summary>
  /// Incomplete factorization LU With Threshold.
  /// </summary>
   static class ILUT
   {
     /// <summary>
     /// Enumerate error values.
     /// </summary>
     typedef enum __eError
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
       } eError;
     
     /// <summary>
     /// Get the error message from the routine <see cref="ILUT"/>.
     /// </summary>
     /// <param name="error">The parameter 'ierr' from the routine.</param>
     /// <returns>The error message related to the error.</returns>
     static const char* ErrorMessage(eError error)
     {
       switch (error)
	 {
	 case Error.Completed:
	   {
	     return
	       "Indicates that the task completed normally.";
	   }
	   
	 case Error.EncounteredAnEmptyRow:
	   {
	     return
	       "Indicates that the routine was interrupted because of an error: the number of elements"
	       + "in some matrix row specified in the sparse format is equal to or less than 0.";
	   }
	   
	 case Error.InvalidDiagonalValue:
	   {
	     return "Indicates that the routine was interrupted because the value of the computed diagonal"
	       +
	       "element is less than the product of the given tolerance and the current matrix row norm, and it cannot be replaced as ipar[30]=0.";
	   }
	   
	 case Error.EncounteredAnInvalidOffset:
	   {
	     return
	       "Indicates that the routine was interrupted because the element ia[i] is less than or equal to the element ia[i - 1] (see Sparse Matrix Storage Format)";
	   }
	   
	 case Error.NotEnoughInternalMemory:
	   {
	     return
	       "Indicates that the routine was interrupted because the memory is insufficient for the internal work arrays.";
	   }
	   
	   
	 case Error.InvalidMaxfil:
	   {
	     return
	       "Indicates that the routine was interrupted because the input value of maxfil is less than 0.";
	   }
		    
		    
	 case Error.InvalidMatrixDimension:
	   {
	     return
	       "Indicates that the routine was interrupted because the size n of the input matrix is less than 0.";
	   }


	 case Error.InvalidMatrixIndexing:
	   {
	     return
	       "Indicates that the routine was interrupted because an element of the array ja is less than 1, or greater than n (see Sparse Matrix Storage Format).";
	   }


	 case Error.InvalidMaxfil2:
	   {
	     return
	       "The value of maxfil is greater than or equal to n. The calculation is performed with the value of maxfil set to (n-1).";
	   }

	 case Error.InvalidTolerance:
	   {
	     return
	       "The value of tol is less than 0. The calculation is performed with the value of the parameter set to (-tol)";
	   }

	 case Error.IncompatibleToleranceAndDiagonalElementThreshold:
	   {
	     return
	       "The absolute value of tol is greater than value of dpar[30]; it can result in instability of the calculation.";
	   }


	 case Error.InvalidDiagonalElementThreshold:
	   {
	     return "The value of dpar[30] is equal to 0. It can cause calculations to fail.";
	   }


	 default:
	   {
	     return "Unknown error value";
	   }
	 }
     };
     
     /// <summary>
     /// Get the required of number of coefficients of the incomplete LU factorization.
     /// </summary>
     /// <param name="n">The dimension.</param>
     /// <param name="maxfil">The maximum fill.</param>
     /// <returns>The required number of coefficients.</returns>
     static I RequiredNumberOfCoefficients(I n,
					   int maxfil)
     {
       return n * (2 * maxfil + 1) - maxfil * (maxfil + 1) + 1;
     };
     
     /// <summary>
     /// ILUT preconditioner based on the incomplete LU factorization with a threshold of a sparse matrix.
     /// </summary>
     /// <remarks>
     /// The routine dcsrilut computes a preconditioner B of a given sparse matrix A stored in the format accepted in the direct sparse solvers:
     /// A~B=L*U , where L is a lower triangular matrix with unit diagonal and U is an upper triangular matrix with non-unit diagonal.
     /// The following threshold criteria are used to generate the incomplete factors L and U:
     /// 1) the resulting entry must be greater than the matrix current row norm multiplied by the parameter <paramref name="tol"/>, and
     /// 2) the number of the non-zero elements in each row of the resulting L and U factors must not be greater than the value of the parameter 
     /// <paramref name="maxfil"/>.
     /// 
     /// </remarks>
     /// 
     /// <param name="n">Size (number of rows or columns) of the original square <paramref name="n"/>-by-<paramref name="n"/> matrix A.</param>
     /// 
     /// <param name="a">Array containing all non-zero elements of the matrix A. The length of the array is equal to their number.</param>
     /// 
     /// <param name="ia">Array of size (<paramref name="n"/>+1) containing indices of non-zero elements in the array <paramref name="a"/>. 
     /// <paramref name="ia"/>[i] is the index of the first non-zero element from the row i. 
     /// The value of the last element <paramref name="ia"/>[n] is equal to the number of non-zeros in the matrix A, plus one.</param>
     /// 
     /// <param name="ja">Array of size equal to the size of the array <paramref name="a"/>. 
     /// This array contains the column numbers for each non-zero element of the matrix A. 
     /// It is important that the indices are in increasing order per row. </param>
     /// 
     /// <param name="bilut">Array containing non-zero elements of the resulting preconditioning matrix B, stored in the format accepted in the direct sparse solvers.
     /// The size of the array is equal to (2*<paramref name="maxfil"/>+1)*<paramref name="n"/>-<paramref name="maxfil"/>*(<paramref name="maxfil"/>+1)+1.
     /// </param>
     /// 
     /// <param name="ibilut">Array of size (n+1) containing indices of non-zero elements in the array bilut. 
     /// <paramref name="ibilut"/>[i] is the index of the first non-zero element from the row i. The value of the last element 
     /// <paramref name="ibilut"/>[<paramref name="n"/>] is equal to the number of non-zeros in the matrix B, plus one.</param>
     /// 
     /// <param name="jbilut">Array, its size is equal to the size of the array <paramref name="bilut"/>. This array contains the column numbers for each non-zero element of the matrix B.</param>
     /// <param name="tol">Tolerance for threshold criterion for the resulting entries of the preconditioner. </param>
     /// <param name="maxfil">Maximum fill-in, which is half of the preconditioner bandwidth. 
     /// The number of non-zero elements in the rows of the preconditioner cannot exceed (2*<paramref name="maxfil"/>+1)</param>
     /// <param name="ipar">Array of size 128. This parameter is used to specify the integer set of data for both the ILUT and RCI FGMRES computations.</param>
     /// <param name="dpar">Array of size 128. This parameter specifies the double precision set of data for both ILUT and RCI FGMRES computations.</param>
     /// <param name="ierr">Error flag, gives information about the routine completion.</param>
     static void Compute(Int32 n,
			 IntPtr a,
			 IntPtr ia,
			 IntPtr ja,
			 void** bilut,
			 void** ibilut,
			 void** jbilut,
			 double tol,
			 Int32 maxfil,
			 Int32[] ipar,
			 double[] dpar,
			 Error* ierr)
     {
       int lierr;
       dcsrilut(&n,
		a,
		ia,
		ja,
		bilut,
		ibilut,
		jbilut,
		&tol,
		&maxfil,
		ipar,
		dpar,
		&lierr);
       
       *ierr = (Error)lierr;
     };
     
     
   };


  
  /// <summary>
  /// Direct solvers.
  /// </summary>
  static class DirectSolvers
  {
    /// <summary>
    /// Parallel Direct Interface Solver.
    /// </summary>
    static class PARDISO
    {
      /// <summary>
      /// Parameters for PARDISO.
      /// </summary>
      class Parameters
      {
      private:
	/// <summary>
	/// Integer parameters.
	/// </summary>
        Int32[] m_iparm = new Int32[64];
	
	/// <summary>
	/// Array of pointers.
	/// </summary>
        IntPtr[] m_pt = new IntPtr[64];
	
	/// <summary>
	/// Get/Set the Fortran indexing.
	/// </summary>
      public:
	bool FortranIndexing
	{
	  get { return m_iparm[34] == 0; }
	  set { m_iparm[34] = (value) ? 0 : 1; }
	}
	
	/// <summary>
	/// Get/Set solve transpose.
	/// </summary>
	bool Transpose
	{
	  get { return m_iparm[11] == 2; }
	  set { m_iparm[11] = (value) ? 2 : 0; }
	}
	
	/// <summary>
	/// The handle of PARDISO.
	/// </summary>
	IntPtr[] Handle
	{
	  get { return m_pt; }
	}
	
	/// <summary>
	/// The arrya of the integer parameters.
	/// </summary>
	Int32[] ParamIntegers
	{
	  get { return m_iparm; }
	}
	
	/// <summary>
	/// Constructor.
	/// </summary>
	Parameters()
	{
	  //
	  // Set the handle to zero.
	  //
	  for (Int32 i = 0; i < m_pt.Length; i++)
	    {
	      m_pt[i] = IntPtr.Zero;
	    }

                 

	  m_iparm[0] = 1; // No solver default
	  m_iparm[1] = 2; // Fill-in reordering from METIS 
	  // Numbers of processors, value of OMP_NUM_THREADS 
	  m_iparm[2] = 32;
	  m_iparm[3] = 31; // No iterative-direct algorithm 
	  m_iparm[4] = 0; // No user fill-in reducing permutation 
	  m_iparm[5] = 0; // Write solution into x 

	  //
	  // OUTPUT:NUmber of iterative refinement steps performed
	  //
	  m_iparm[6] = 0; // Not in use 


	  m_iparm[7] = 2; // Max numbers of iterative refinement steps 
	  m_iparm[8] = 0; // Not in use 
	  m_iparm[9] = 13; // Perturb the pivot elements with 1E-13 

	  //
	  // Scaling factors
	  //
	  m_iparm[10] = 1; // Use nonsymmetric permutation and scaling MPS 

	  m_iparm[12] = 1;
	  // Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy 

	  m_iparm[13] = 0; // Output: Number of perturbed pivots 
	  m_iparm[14] = 0; // Not in use 
	  m_iparm[15] = 0; // Not in use 
	  m_iparm[16] = 0; // Not in use 
	  m_iparm[17] = -1; // Output: Number of nonzeros in the factor LU 
	  m_iparm[18] = -1; // Output: Mflops for LU factorization 
	  m_iparm[19] = 0; // Output: Numbers of CG Iterations 



	  //
	  // We choose C-indexing by default
	  //
	  FortranIndexing = false;

	  //
	  // We solve Ax=b by default
	  //
	  Transpose = false;

	  //
	  // Reordering
	  //
	  Reordering = MatrixReordering.MinimumDegree;

	}
	
	/// <summary>
	/// Enumerate options for the matrix reordering.
	/// </summary>
	enum MatrixReordering
	  {
	    /// <summary>
	    /// 
	    /// </summary>
	    MinimumDegree = 0,
	    /// <summary>
	    /// 
	    /// </summary>
	    NestedDissectionFromMetis = 2,
	    /// <summary>
	    /// 
	    /// </summary>
	    OpenMPVersionOfNestedDissection = 3
	  };
	
	
	/// <summary>
	/// Get/Set the matrix reordering to apply.
	/// </summary>
	MatrixReordering Reordering
	{
	  get { return (MatrixReordering)m_iparm[1]; }
	  
	  set
	    {
	      ValidationUtil.ValidateEnumArg(value, typeof(MatrixReordering), "value");
	      m_iparm[1] = (int)value;
	    }
	};
	
      };
	
      /// <summary>
      /// Enumerate error values.
      /// </summary>
      enum Error
      {
               /// <summary>
               /// no error.
               /// </summary>
               NoError = 0,

               /// <summary>
               /// input inconsistency.
               /// </summary>
               InconsistentInput = -1,

               /// <summary>
               /// not enough memory.
               /// </summary>
               NotEnoughMemory = -2,
               
               /// <summary>
               /// reordering problem.
               /// </summary>
               ReorderingProblem = -3,

               /// <summary>
               /// zero pivot, numerical factorization or iterative refinement problem.
               /// </summary>
               ZeroPivotNumericalFactorizationOrRefinementProblem = -4,

               /// <summary>
               /// unclassified (internal) error.
               /// </summary>
               Unclassified = -5,

               /// <summary>
               /// pre-ordering failed (matrix types 11, 13 only).
               /// </summary>
               PreorderingFailed = -6,

               /// <summary>
               /// diagonal matrix is singular.
               /// </summary>
               SingularDiagonalMatrix = -7,

               /// <summary>
               /// 32-bit integer overflow problem.
               /// </summary>
               Int32Overflow = -8,

               /// <summary>
               /// not enough memory for OOC.
               /// </summary>
               NotEnoughMemoryForOOC = -9,
               
               /// <summary>
               /// problems with opening OOC temporary file.
               /// </summary>
               ProblemsWithOpeningOOCTemporaryFiles = -10,

               /// <summary>
               /// read/write problems with the OOC data file.
               /// </summary>
               ReadWriteProblemsWithTheOOCDataFile = -11
            }

            /// <summary>
            /// Get the related error message.
            /// </summary>
            /// <param name="error"></param>
            /// <returns></returns>
            public static string GetErrorMessage(Error error)
            {
               ValidationUtil.ValidateEnumArg(error, typeof(Error), "error");
               switch (error)
               {
                  case Error.InconsistentInput:
                     {
                        return "input inconsistent";
                     }
                  case Error.NotEnoughMemory:
                     {
                        return "not enough memory";
                     }
                  case Error.ReorderingProblem:
                     {
                        return "reordering problem";
                     }
                  case Error.ZeroPivotNumericalFactorizationOrRefinementProblem:
                     {
                        return "zero pivot, numerical factorization or iterative refinement problem";
                     }
                  case Error.Unclassified:
                     {
                        return "unclassified (internal) error";
                     }
                  case Error.PreorderingFailed:
                     {
                        return "pre-ordering failed (matrix types 11, 13 only)";
                     }
                  case Error.SingularDiagonalMatrix:
                     {
                        return "diagonal matrix is singular";
                     }
                  case Error.Int32Overflow:
                     {
                        return "32-bit integer overflow problem";
                     }
                  case Error.NotEnoughMemoryForOOC:
                     {
                        return "not enough memory for OOC";
                     }
                  case Error.ProblemsWithOpeningOOCTemporaryFiles:
                     {
                        return "problems with opening OOC temporary file";
                     }
                  case Error.ReadWriteProblemsWithTheOOCDataFile:
                     {
                        return "read/write problems with the OOC data file";
                     }

                  default:
                     {
                        Debug.Assert(Error.NoError == error, "Unexpected value of the enumeration Error.");
                        return "no error";
                     }
               }
            }

            /// <summary>
            /// Enumerates job values.
            /// </summary>
            public enum Job
            {
               /// <summary>
               /// Release the memory.
               /// </summary>
               ReleaseMemory = -1,

               /// <summary>
               /// Perform symbolic factorization.
               /// </summary>
               SymbolicFactorization = 11,
               
               /// <summary>
               /// Perform numerical factorization.
               /// </summary>
               NumericalFactorization = 22,

               /// <summary>
               /// Solve.
               /// </summary>
               Solve = 33,
            }


            /// <summary>
            /// Enumerate the matrix types.
            /// </summary>
            public enum MatrixType
            {
               /// <summary>
               /// Real and structurally symmetric.
               /// </summary>
               RealAndStructurallySymmetric = 1,

               /// <summary>
               /// Real and symmetric positive definite.
               /// </summary>
               RealAndSymmetricPositiveDefinite = 2,

               /// <summary>
               /// Real and symmetric indefinite.
               /// </summary>
               RealAndSymmetricIndefinite = -2,

               /// <summary>
               /// Real and not symmetric.
               /// </summary>
               ReadAndUnsymmetric = 11
            }

            /// <summary>
            /// Routine to call the PARDISO solver.
            /// </summary>
            /// <param name="parameters">The parameters.</param>
            /// <param name="maxfct">Maximal number of factors in memory.</param>
            /// <param name="mnum">The number of matrix (from 1 to <paramref name="maxfct"/>) to solve.</param>
            /// <param name="mtype">Matrix type.</param>
            /// <param name="job">Controls the execution of the solver.</param>
            /// <param name="n">Number of equations in the sparse linear system A*X = B.</param>
            /// <param name="a">Contains the non-zero elements of the coefficient matrix A.</param>
            /// <param name="ia">rowIndex array in CSR format.</param>
            /// <param name="ja">columns array in CSR format.</param>
            /// <param name="perm">Holds the permutation vector of size <paramref name="n"/>.</param>
            /// <param name="nrhs">Number of right-hand sides that need to be solved for.</param>
            /// <param name="msglvl">Message level information.</param>
            /// <param name="b">Right hand side vectors.</param>
            /// <param name="x">Solution vectors.</param>
            /// <param name="outError">OUTPUT:Error indicator.</param>
            public static unsafe void Pardiso(Parameters parameters,
                                              Int32 maxfct,
                                              Int32 mnum,
                                              Int32 mtype,
                                              Job job,
                                              Int32 n,
                                              [In] IntPtr a,
                                              [In] IntPtr ia,
                                              [In] IntPtr ja,
                                              [In] IntPtr perm,
                                              [In] Int32 nrhs,
                                              [In] Int32 msglvl,
                                              [In, Out] IntPtr b,
                                              [Out] IntPtr x,
                                              [Out] out Error outError)
            {
               ValidationUtil.ValidateEnumArg(job, typeof(Job), "job");

               int err;
               int phase = (int) job;

               pardiso(parameters.Handle,
                       &maxfct,
                       &mnum,
                       &mtype,
                       &phase,
                       &n,
                       a,
                       ia,
                       ja,
                       perm,
                       &nrhs,
                       parameters.ParamIntegers,
                       &msglvl,
                       b,
                       x,
                       out err);

               outError = (Error) err;
            }


            [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true,
               SetLastError = false)]
            private static extern unsafe void pardiso([In, Out] IntPtr[] handle,
                                                      [In] Int32* maxfct,
                                                      [In] Int32* mnum,
                                                      [In] Int32* mtype,
                                                      Int32* phase,
                                                      [In] Int32* n,
                                                      [In] IntPtr a,
                                                      [In] IntPtr ia,
                                                      [In] IntPtr ja,
                                                      [In] IntPtr perm,
                                                      [In] Int32* nrhs,
                                                      [In, Out] Int32[] iparm,
                                                      [In] Int32* msglvl,
                                                      [In, Out] IntPtr b,
                                                      [Out] IntPtr x,
                                                      [Out] out Int32 error);
         }

      }

      #endregion

      #region Iterative solvers

      /// <summary>
      /// Iterative sparse solvers.
      /// </summary>
      public static class IterativeSparseSolvers
      {
         /// <summary>
         /// Parameters of the Conjugate Gradient Solver.
         /// </summary>
         public class ParametersBase
         {

            /// <summary>
            /// The array of integer parameters.
            /// </summary>
            protected readonly Int32[] m_ipar = new Int32[128];

            /// <summary>
            /// The array of double parameters.
            /// </summary>
            protected readonly double[] m_dpar = new double[128];

            /// <summary>
            /// The array of the integer parameters.
            /// </summary>
            public Int32[] ParamIntegers
            {
               get { return m_ipar; }
            }

            /// <summary>
            /// The array of the real parameters.
            /// </summary>
            public double[] ParamReals
            {
               get { return m_dpar; }
            }


         }


         /// <summary>
         /// Flexible GMRES
         /// </summary>
         public static class FGMRES
         {

            /// <summary>
            /// Parameters of the MKL FGMRES Solver.
            /// </summary>
            public class Parameters : ParametersBase
            {
               /// <summary>
               /// Print information.
               /// </summary>
               /// <param name="textWriter">The text writer.</param>
               [CoverageExclude("Categorized with Conditional(\"DEBUG\")")]
               [Conditional("DEBUG")]
               public void PrintInformation(TextWriter textWriter)
               {
                  ValidationUtil.ValidateArgNotNull(textWriter, "textWriter");

                  textWriter.WriteLine("Some info about the current run of RCI FGMRES method:");
                  if (m_ipar[7] != 0)
                  {
                     textWriter.WriteLine(
                        "   As ipar[7]={0}, the automatic test for the maximal number of iterations will be", m_ipar[7]);
                     textWriter.WriteLine("   performed");
                  }
                  else
                  {
                     textWriter.WriteLine(
                        "   As ipar[7]={0}, the automatic test for the maximal number of iterations will be", m_ipar[7]);
                     textWriter.WriteLine("   skipped");
                  }
                  textWriter.WriteLine("+++");
                  if (0 != m_ipar[8])
                  {
                     textWriter.WriteLine("   As ipar[8]={0}, the automatic residual test will be performed", m_ipar[8]);
                  }
                  else
                  {
                     textWriter.WriteLine("   As ipar[8]={0}, the automatic residual test will be skipped", m_ipar[8]);
                  }
                  textWriter.WriteLine("+++");
                  if (0 != m_ipar[9])
                  {
                     textWriter.WriteLine("   As ipar[9]={0}, the user-defined stopping test will be requested via",
                                          m_ipar[9]);
                     textWriter.WriteLine("   RCI_request=2");
                  }
                  else
                  {
                     textWriter.WriteLine(
                        "   As ipar[9]={0}, the user-defined stopping test will not be requested, thus,", m_ipar[9]);
                     textWriter.WriteLine("   RCI_request will not take the value 2");
                  }
                  textWriter.WriteLine("+++");
                  if (0 != m_ipar[10])
                  {
                     textWriter.WriteLine(
                        "   As ipar[10]={0}, the Preconditioned FGMRES iterations will be performed, thus,", m_ipar[10]);
                     textWriter.WriteLine("   the preconditioner action will be requested via RCI_request=3");
                  }
                  else
                  {
                     textWriter.WriteLine(
                        "   As ipar[10]={0}, the Preconditioned FGMRES iterations will not be performed,", m_ipar[10]);
                     textWriter.WriteLine("   thus, RCI_request will not take the value 3");
                  }
                  textWriter.WriteLine("+++");
                  if (0 != m_ipar[11])
                  {
                     textWriter.WriteLine(
                        "   As ipar[11]={0}, the automatic test for the norm of the next generated vector is",
                        m_ipar[11]);
                     textWriter.WriteLine(
                        "   not equal to zero up to rounding and computational errors will be performed,");
                     textWriter.WriteLine("   thus, RCI_request will not take the value 4");
                  }
                  else
                  {
                     textWriter.WriteLine(
                        "   As ipar[11]={0}, the automatic test for the norm of the next generated vector is",
                        m_ipar[11]);
                     textWriter.WriteLine(
                        "   not equal to zero up to rounding and computational errors will be skipped,");
                     textWriter.WriteLine("   thus, the user-defined test will be requested via RCI_request=4");
                  }
                  textWriter.WriteLine("+++");
               }


               /// <summary>
               /// Set the relative tolerance.
               /// </summary>
               public double RelativeTolerance
               {
                  set { m_dpar[0] = value; }
               }

               /// <summary>
               /// Make the residual stopping test automatic.
               /// </summary>
               public bool AutomaticResidualStoppingTest
               {
                  set { m_ipar[8] = value ? 1 : 0; }
               }

               /// <summary>
               /// Set the maximum number of iterations.
               /// </summary>
               public Int32 MaximumNumberOfIterations
               {
                  set { m_ipar[4] = value; }
               }

               /// <summary>
               /// Set the maximum number of non-restarted iterations.
               /// </summary>
               public Int32 MaximumNumberOfNonRestartedIterations
               {
                  set { m_ipar[14] = value; }
               }

               /// <summary>
               /// Get the norm of the orthogonal vector.
               /// </summary>
               public double NormOfOrthogonalVector
               {
                  get { return m_dpar[6]; }
               }

               /// <summary>
               /// Get/Set the threshold of the norm of the orthogonal vector.
               /// </summary>
               public double ZeroNormThresholdOfOrthogonalVector
               {
                  get { return m_dpar[7]; }
                  set { m_dpar[7] = value; }
               }

               /// <summary>
               /// Get the index of the source temporary vector.
               /// </summary>
               public Int32 IndexOfSourceTemporaryVector
               {
                  get { return m_ipar[21] - 1; }
               }

               /// <summary>
               /// Get the index of output temporary vector.
               /// </summary>
               public Int32 IndexOfOutputTemporaryVector
               {
                  get { return m_ipar[22] - 1; }
               }

               /// <summary>
               /// Make the user defined stopping test automatic.
               /// </summary>
               public bool AutomaticUserDefinedStoppingTest
               {
                  set { m_ipar[9] = value ? 1 : 0; }
               }

               /// <summary>
               /// Make the zero norm of the orthogonal vector automatic.
               /// </summary>
               public bool AutomaticTestForZeroNormOfOrthogonalVector
               {
                  get { return m_ipar[11] != 0; }
                  set { m_ipar[11] = value ? 1 : 0; }
               }

               /// <summary>
               /// Indicates if the method is preconditioned.
               /// </summary>
               public bool PreconditionedMethod
               {
                  get { return m_ipar[10] == 1; }
                  set { m_ipar[10] = value ? 1 : 0; }
               }

               /// <summary>
               /// The required memory.
               /// </summary>
               /// <param name="maximumNumberOfNonRestartedIterations">The maximum number of non-restarted iterations.</param>
               /// <param name="n">The dimension.</param>
               /// <returns>The required memory.</returns>
               public static Int32 GetRequiredMemory(Int32 maximumNumberOfNonRestartedIterations, Int32 n)
               {
                  return ((2 * maximumNumberOfNonRestartedIterations + 1) * n +
                          maximumNumberOfNonRestartedIterations * (maximumNumberOfNonRestartedIterations + 9) / 2 + 1);
               }
            }

            /// <summary>
            /// Enumerate the values of the parameter request.
            /// </summary>
            public enum Request
            {
               /// <summary>
               /// Indicates that the task completed normally and the solution is found and stored in the vector x. 
               /// This occurs only if the stopping tests are fully automatic. 
               /// For the user defined stopping tests, see the description of the RCI_request= 2 or 4.
               /// </summary>
               Completed = 0,

               /// <summary>
               /// Indicates that the routine was interrupted because the maximum number of iterations was reached, but the relative stopping criterion was not met. 
               /// This situation occurs only if you request both tests.
               /// </summary>
               MaximumNumberOfIterationsReached = -1,

               /// <summary>
               /// Indicates that the routine was interrupted because of an attempt to divide by zero. 
               /// Usually this happens if the matrix is degenerate or almost degenerate. 
               /// However, it may happen if the parameter dpar is altered, or if the method is not stopped when the solution is found.
               /// </summary>
               AttemptDivideByZero = -10,

               /// <summary>
               /// Indicates that the routine was interrupted because it entered an infinite cycle. 
               /// Usually this happens because the values ipar[7], ipar[8], ipar[9] were altered outside of the routine, 
               /// or the dfgmres_check routine was not called.
               /// </summary>
               InfiniteCycle = -11,

               /// <summary>
               /// Indicates that the routine was interrupted because errors were found in the method parameters. 
               /// Usually this happens if the parameters ipar and dpar were altered by mistake outside the routine.
               /// </summary>
               InconsitentParameters = -12,

               /// <summary>
               /// Indicates that you must multiply the matrix by rwork[ipar[21] - 1:ipar[21] + n - 2], 
               /// put the result in the rwork[ipar[22] - 1:ipar[22] + n - 2], and return control back to the routine dfgmres.
               /// </summary>
               ApplyMatrixVectorProduct = 1,

               /// <summary>
               /// Indicates that you must perform the stopping tests. If they fail, return control to the dfgmres routine. Otherwise, the FGMRES solution is found, and you can run the fgmres_get routine to update the computed solution in the vector x.
               /// </summary>
               ApplyStoppingTests = 2,

               /// <summary>
               /// Indicates that you must apply the inverse preconditioner to tmp[ipar[21] - 1:ipar[21] + n - 2], put the result in the tmp[ipar[22] - 1:ipar[22] + n - 2], and return control back to the routine dfgmres.
               /// </summary>
               ApplyPreconditioner = 3,

               /// <summary>
               /// Indicates that you must check the norm of the currently generated vector. 
               /// If it is not zero within the computational/rounding errors, return control to the dfgmres routine. 
               /// Otherwise, the FGMRES solution is found, and you can run the dfgmres_get routine to update the computed solution in the vector x.
               /// </summary>
               CheckTheNormOfTheGeneratedVector = 4
            }

            /// <summary>
            /// Get the error message related to the status.
            /// </summary>
            /// <param name="request">The parameter 'request' from the routine.</param>
            public static string ErrorMessage(Request request)
            {
               ValidationUtil.ValidateEnumArg(request, typeof(Request), "request");
               switch (request)
               {
                  case Request.MaximumNumberOfIterationsReached:
                     {
                        return
                           "Indicates that the routine was interrupted because the maximum number of iterations was reached, but the relative stopping criterion was not met."
                           + "This situation occurs only if you request both tests.";
                     }

                  case Request.AttemptDivideByZero:
                     {
                        return
                           "Indicates that the routine was interrupted because of an attempt to divide by zero. "
                           + "Usually this happens if the matrix is degenerate or almost degenerate. "
                           +
                           "However, it may happen if the parameter dpar is altered, or if the method is not stopped when the solution is found.";
                     }

                  case Request.InfiniteCycle:
                     {
                        return "Indicates that the routine was interrupted because it entered an infinite cycle. "
                               +
                               "Usually this happens because the values ipar[7], ipar[8], ipar[9] were altered outside of the routine,"
                               + "or the dfgmres_check routine was not called.";
                     }

                  case Request.InconsitentParameters:
                     {
                        return "Indicates that the routine was interrupted because errors were found in the method parameters. "
                               +
                               "Usually this happens if the parameters ipar and dpar were altered by mistake outside the routine.";
                     }

                  case Request.Completed:
                     {
                        return "Completed";
                     }
                  case Request.ApplyMatrixVectorProduct:
                     {
                        return "ApplyMatrixVectorProduct";
                     }
                  case Request.ApplyStoppingTests:
                     {
                        return "ApplyStoppingTests";
                     }
                  case Request.ApplyPreconditioner:
                     {
                        return "ApplyPreconditioner";
                     }
                  case Request.CheckTheNormOfTheGeneratedVector:
                     {
                        return "ApplyPreconditioner";
                     }
                  default:
                     {
                        int irequest = (int) request;
                        if (irequest < 0)
                        {
                           switch (irequest)
                           {
                              case -10000:
                                 {
                                    return "Indicates that the routine failed to complete the task.";
                                 }
                              default:
                                 {
                                    return "Unknown error";
                                 }
                           }
                        }
                        else
                        {
                           return "Unknown request value";
                        }
                     }
               }
            }

            #region Init

            /// <summary>
            /// Initializes the solver.
            /// </summary>
            /// <param name="n">INPUT:Sets the size of the problem.</param>
            /// <param name="request">OUTPUT:Gives information about the result of the routine.</param>
            /// <param name="parameters">OUTPUT: The parameters. </param>
            /// <param name="rwork">OUTPUT: Array of size ((2*<paramref name="parameters"/>[14] + 1)*<paramref name="n"/> + <paramref name="parameters"/>[14]*(<paramref name="parameters"/>[14 - 1] + 9)/2 + 1).</param>
            public static void Init([In] Int32 n,
                                    [In, Out] ref Request request,
                                    Parameters parameters,
                                    [Out] IntPtr rwork)
            {
               Int32 ln = n;
               Int32 lrequest = (Int32) request;

               dfgmres_init(n : ref ln,
                            sol : IntPtr.Zero,
                            rhs : IntPtr.Zero,
                            request : ref lrequest,
                            ipar : parameters.ParamIntegers,
                            dpar : parameters.ParamReals,
                            rwork : rwork);

               request = (Request) lrequest;
            }

            [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true,
               SetLastError = false)]
            private static extern void dfgmres_init([In, Out] ref Int32 n,
                                                    [In] IntPtr sol,
                                                    [In] IntPtr rhs,
                                                    [In, Out] ref Int32 request,
                                                    [Out] Int32[] ipar,
                                                    [Out] double[] dpar,
                                                    [Out] IntPtr rwork);

            #endregion

            #region Check

            /// <summary>
            /// Checks consistency and correctness of the user defined data.
            /// </summary>
            /// <param name="parameters">OUTPUT: The parameters. </param>
            /// <param name="n">INPUT:Sets the size of the problem.</param>
            /// <param name="rwork">OUTPUT: Array of size ((2*<paramref name="parameters"/>[14] + 1)*<paramref name="n"/> + <paramref name="parameters"/>[14]*(<paramref name="parameters"/>[14] + 9)/2 + 1).</param>       
            /// <param name="request">OUTPUT: Gives information about result of the routine.</param>
            public static void Check(Parameters parameters,
                                     [In] Int32 n,
                                     [Out] IntPtr rwork,
                                     [In, Out] ref Request request)
            {
               Int32 ln = n;
               Int32 lrequest = (Int32) request;

               dfgmres_check(n : ref ln,
                             sol : IntPtr.Zero,
                             rhs : IntPtr.Zero,
                             request : ref lrequest,
                             ipar : parameters.ParamIntegers,
                             dpar : parameters.ParamReals,
                             rwork : rwork);

               request = (Request) lrequest;
            }


            [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true,
               SetLastError = false)]
            private static extern void dfgmres_check([In, Out] ref Int32 n,
                                                     [In] IntPtr sol,
                                                     [In] IntPtr rhs,
                                                     [In, Out] ref Int32 request,
                                                     [Out] Int32[] ipar,
                                                     [Out] double[] dpar,
                                                     [Out] IntPtr rwork);

            #endregion

            #region Get

            /// <summary>
            /// Retrieves the number of the current iteration and updates the solution.
            /// </summary>
            /// <param name="parameters">OUTPUT: The parameters. </param>
            /// <param name="n">INPUT:Sets the size of the problem.</param>
            /// <param name="sol">OUTPUT: Array of size <paramref name="n"/>. Contains the initial approximation to the solution vector. Normally it is equal to 0 or to <paramref name="rhs"/>.</param>
            /// <param name="rhs">OUTPUT: Array of size <paramref name="n"/>. Contains the right-hand side vector.</param>
            /// <param name="request">OUTPUT: Gives information about result of the routine.</param>
            /// <param name="rwork">INPUT:Array of size ((2*<paramref name="parameters"/>[14] + 1)*<paramref name="n"/> + <paramref name="parameters"/>[14]*(<paramref name="parameters"/>[14] + 9)/2 + 1).</param>
            /// <param name="itercount">OUTPUT: Contains the value of the current iteration number.</param>
            public static  void Get(Parameters parameters,
                                          [In] Int32 n,
                                          [Out] IntPtr sol,
                                          [Out] IntPtr rhs,
                                          [In] IntPtr rwork,
                                          [In, Out] ref Request request,
                                          [In, Out] ref Int32 itercount)
            {
               Int32 ln = n;
               Int32 lrequest = (Int32) request;

               dfgmres_get(ref ln,
                           sol,
                           rhs,
                           ref lrequest,
                           parameters.ParamIntegers,
                           parameters.ParamReals,
                           rwork,
                           ref itercount);

               request = (Request) lrequest;

            }

            [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true,
               SetLastError = false)]
            private static extern void dfgmres_get([In, Out] ref Int32 n,
                                                   [Out] IntPtr sol,
                                                   [Out] IntPtr rhs,
                                                   [In, Out] ref Int32 request,
                                                   [In] Int32[] ipar,
                                                   [In] double[] dpar,
                                                   [In] IntPtr rwork,
                                                   [In, Out] ref Int32 itercount);

            #endregion

            #region Run

            /// <summary>
            /// Makes the FGMRES iterations.
            /// </summary>
            /// <param name="parameters">OUTPUT: The parameters. </param>
            /// <param name="n">INPUT:Sets the size of the problem.</param>
            /// <param name="sol">INPUT: Array of size <paramref name="n"/>. Contains the initial approximation to the solution vector. Normally it is equal to 0 or to b.</param>
            /// <param name="rhs">INPUT: Array of size <paramref name="n"/>. Contains the right-hand side vector.</param>
            /// <param name="rwork">INPUT/OUTPUT:The working array.</param>
            /// <param name="request">OUTPUT: Gives information about result of the routine.</param>
            public static void Run(Parameters parameters,
                                          [In] Int32 n,
                                          [In] IntPtr sol,
                                          [In] IntPtr rhs,
                                          [In, Out] IntPtr rwork,
                                          [In, Out] ref Request request)
            {
               Debug.Assert(null != parameters, "The parameter 'parameters' cannot be null.");
               Debug.Assert(null != rwork, "The parameter 'rwork' cannot be null.");
               Debug.Assert(null != sol, "The parameter 'sol' cannot be null.");
               Debug.Assert(null != rhs, "The parameter 'rhs' cannot be null.");

               Int32 ln = n;
               Int32 lrequest = (Int32) request;
               dfgmres(ref ln,
                       sol,
                       rhs,
                       ref lrequest,
                       parameters.ParamIntegers,
                       parameters.ParamReals,
                       rwork);
               request = (Request) lrequest;
            }


            [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true,
               SetLastError = false)]
            private static extern void dfgmres(ref Int32 n,
                                                      IntPtr sol,
                                                      IntPtr rhs,
                                                      ref Int32 request,
                                                      Int32[] ipar,
                                                      double[] dpar,
                                                      IntPtr rwork);

            #endregion

         }

         /// <summary>
         /// Conjugate gradient.
         /// </summary>
         public static class CG
         {

            /// <summary>
            /// Parameters of the Conjugate Gradient Solver.
            /// </summary>
            public class Parameters : ParametersBase
            {

               /// <summary>
               /// If the value is not equal to 0, the dcg/dcgmrhs routine performs the residual stopping test: dpar(5) ≤ dpar(4)= dpar(1)*dpar(3)+ dpar(2). 
               /// Otherwise, the method is stopped and corresponding value is assigned to the RCI_request.
               /// If the value is 0, the routine does not perform this stopping test. The default value is 0.
               /// </summary>
               public bool AutomaticResidualStoppingTest
               {
                  get { return m_ipar[8] == 1; }
                  set { m_ipar[8] = value ? 1 : 0; }
               }

               /// <summary>
               /// Set the maximum number of iterations.
               /// </summary>
               public Int32 MaximumNumberOfIterations
               {
                  set { m_ipar[4] = value; }
               }


               /// <summary>
               /// If the value is not equal to 0, the dcg/dcgmrhs routine 
               /// requests a user-defined stopping test by setting the output parameter RCI_request=2. 
               /// If the value is 0, the routine does not perform the user defined stopping test. The default value is 1.
               /// </summary>
               public bool AutomaticUserDefinedStoppingTest
               {
                  get { return m_ipar[9] == 1; }
                  set { m_ipar[9] = value ? 1 : 0; }
               }

               /// <summary>
               /// The value is equal to 0, the dcg/dcgmrhs routine runs the non-preconditioned version of the corresponding CG method. 
               /// Otherwise, the routine runs the preconditioned version of the CG method, and by setting the output parameter RCI_request=3, 
               /// indicates that you must perform the preconditioning step. The default value is 0.
               /// </summary>
               public bool PreconditionedMethod
               {
                  get { return m_ipar[10] == 1; }
                  set { m_ipar[10] = value ? 1 : 0; }
               }

               /// <summary>
               /// Set the relative tolerance.
               /// </summary>
               /// <remarks>Default value is 1e-6.</remarks>
               public double RelativeTolerance
               {
                  set { m_dpar[0] = value; }
               }

               /// <summary>
               /// Set the absolute tolerance.
               /// </summary>
               /// <remarks>Default value is 0.0.</remarks>
               public double AbsoluteTolerance
               {
                  set { m_dpar[1] = value; }
               }

               /// <summary>
               /// The required memory.
               /// </summary>
               /// <param name="n">The dimension.</param>
               /// <returns>The required memory.</returns>
               public static Int32 GetRequiredMemory(Int32 n)
               {
                  return 4 * n;
               }
            }

            /// <summary>
            /// Enumerate the values of the parameter request.
            /// </summary>
            public enum Request
            {
               /// <summary>
               /// Indicates that the task completed normally and the solution is found and stored in the vector x. 
               /// This occurs only if the stopping tests are fully automatic. 
               /// For the user defined stopping tests, see the description of the RCI_request= 2 or 4.
               /// </summary>
               Completed = 0,

               /// <summary>
               /// Indicates that the routine was interrupted because the maximum number of iterations was reached, but the relative stopping criterion was not met. 
               /// This situation occurs only if you request both tests.
               /// </summary>
               RunMaximumNumberOfIterationsReached = -1,

               /// <summary>
               /// Indicates that the routine was interrupted because of an attempt to divide by zero. 
               /// This situation happens if the matrix is non-positive definite or almost non-positive definite.
               /// </summary>
               RunAttemptDivideByZero = -2,

               /// <summary>
               /// Indicates that the routine was interrupted because the residual norm is invalid. 
               /// This usually happens because the value dpar(6) was altered outside of the routine, 
               /// or the dcg_check routine was not called.
               /// </summary>
               RunInvalidResidualNorm = -10,

               /// <summary>
               /// Indicates that the routine was interrupted because it entered an infinite cycle. 
               /// Usually this happens because the values ipar[7], ipar[8], ipar[9] were altered outside of the routine, 
               /// or the dfgmres_check routine was not called.
               /// </summary>
               RunInfiniteCycle = -11,

               /// <summary>
               /// Multiply the matrix by tmp(1:n,1), put the result in tmp(1:n,2), and return the control to the dcg/dcgmrhs routine;
               /// </summary>
               RunApplyMatrixVectorProduct = 1,

               /// <summary>
               /// Indicates that you must perform the stopping tests. 
               /// If they fail, return the control to the dcg/dcgmrhs routine. 
               /// If the stopping tests succeed, it indicates that the solution is found and stored in the x array;
               /// </summary>
               RunApplyStoppingTests = 2,

               /// <summary>
               /// Indicates that you must apply  the preconditioner to tmp(1:n,3), put the result in tmp(1:n,4), and return the control to the dcg routine;
               /// </summary>
               RunApplyPreconditioner = 3,

               /// <summary>
               /// Failure to complete the task.
               /// </summary>
               InitFailure = -10000,

               /// <summary>
               /// Indicates that the task is interrupted and the errors occur.
               /// </summary>
               CheckInterrupted = -1100,

               /// <summary>
               /// Indicates that there are some warning messages.
               /// </summary>
               CheckWarningMessages = -1001,

               /// <summary>
               /// Indicates that the routine changed some parameters to make them consistent or correct.
               /// </summary>
               CheckCorrectedParameters = -1010,

               /// <summary>
               /// Indicates that there are some warning messages and that the routine changed some parameters.
               /// </summary>
               CheckWarningMessagesAndCorrectedParameters = -1011,


            }

            /// <summary>
            /// Get the error message related to the status.
            /// </summary>
            /// <param name="request">The parameter 'request' from the routine.</param>
            public static string ErrorMessage(Request request)
            {
               ValidationUtil.ValidateEnumArg(request, typeof(Request), "request");
               switch (request)
               {
                  case Request.RunMaximumNumberOfIterationsReached:
                     {
                        return
                           "Indicates that the routine was interrupted because the maximum number of iterations was reached, but the relative stopping criterion was not met."
                           + "This situation occurs only if you request both tests.";
                     }

                  case Request.RunAttemptDivideByZero:
                     {
                        return
                           "Indicates that the routine was interrupted because of an attempt to divide by zero. "
                           + "Usually this happens if the matrix is degenerate or almost degenerate. "
                           +
                           "However, it may happen if the parameter dpar is altered, or if the method is not stopped when the solution is found.";
                     }

                  case Request.RunInfiniteCycle:
                     {
                        return "Indicates that the routine was interrupted because it entered an infinite cycle. "
                               +
                               "Usually this happens because the values ipar[7], ipar[8], ipar[9] were altered outside of the routine,"
                               + "or the dfgmres_check routine was not called.";
                     }

                  case Request.InitFailure:
                     {
                        return "Failure to complete the task.";
                     }

                  case Request.CheckInterrupted:
                     {
                        return " Indicates that the task is interrupted and the errors occur.";
                     }

                  case Request.CheckWarningMessages:
                     {
                        return " Indicates that there are some warning messages.";
                     }

                  case Request.CheckCorrectedParameters:
                     {
                        return " Indicates that the routine changed some parameters to make them consistent or correct.";
                     }

                  case Request.CheckWarningMessagesAndCorrectedParameters:
                     {
                        return
                           " Indicates that there are some warning messages and that the routine changed some parameters.";
                     }

                  case Request.Completed:
                     {
                        return "Completed";
                     }
                  case Request.RunApplyMatrixVectorProduct:
                     {
                        return "ApplyMatrixVectorProduct";
                     }
                  case Request.RunApplyStoppingTests:
                     {
                        return "ApplyStoppingTests";
                     }
                  case Request.RunApplyPreconditioner:
                     {
                        return "ApplyPreconditioner";
                     }

                  default:
                     {
                        return "Unknown request value";
                     }
               }
            }

            #region Init

            /// <summary>
            /// Initializes the solver.
            /// </summary>
            /// <param name="n">INPUT:Sets the size of the problem.</param>
            /// <param name="request">OUTPUT:Gives information about the result of the routine.</param>
            /// <param name="parameters">OUTPUT: The parameters. </param>
            /// <param name="rwork">OUTPUT:  Array of size (<paramref name="n"/>,4).</param>
            public static void Init([In] Int32 n,
                                    [In, Out] ref Request request,
                                    Parameters parameters,
                                    [Out] IntPtr rwork)
            {
               Int32 ln = n;
               Int32 lrequest = (Int32) request;
               
               dcg_init(ref ln,
                        sol : IntPtr.Zero, // no need of the solution vector
                        rhs : IntPtr.Zero, // no need of the right-hand side vector
                        request : ref lrequest,
                        ipar : parameters.ParamIntegers,
                        dpar : parameters.ParamReals,
                        rwork : rwork);

               request = (Request) lrequest;
            }


            [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true,
               SetLastError = false)]
            private static extern void dcg_init([In, Out] ref Int32 n,
                                                [In] IntPtr sol,
                                                [In] IntPtr rhs,
                                                [In, Out] ref Int32 request,
                                                [Out] Int32[] ipar,
                                                [Out] double[] dpar,
                                                [Out] IntPtr rwork);

            #endregion

            #region Check

            /// <summary>
            /// Checks consistency and correctness of the user defined data.
            /// </summary>
            /// <param name="parameters">OUTPUT: The parameters. </param>
            /// <param name="n">INPUT:Sets the size of the problem.</param>
            /// <param name="rwork">OUTPUT: Array of size (<paramref name="n"/>,4)).</param>       
            /// <param name="request">OUTPUT: Gives information about result of the routine.</param>
            /// <remarks>
            /// The routine dcg_check checks consistency and correctness of the parameters to be passed to the solver routine dcg. 
            /// However this operation does not guarantee that the solver returns the correct result. It only reduces the chance 
            /// of making a mistake in the parameters of the method. Skip this operation only if you are sure that the correct data 
            /// is specified in the solver parameters. The lengths of all vectors must be defined in a previous call to the dcg_init 
            /// routine.
            ///</remarks>
            public static void Check(Parameters parameters,
                                     [In] Int32 n,
                                     [Out] IntPtr rwork,
                                     [In, Out] ref Request request)
            {
               Int32 ln = n;
               Int32 lrequest = (Int32) request;

               dcg_check(n : ref ln,
                         sol : IntPtr.Zero,
                         rhs: IntPtr.Zero,
                         request : ref lrequest,
                         ipar : parameters.ParamIntegers,
                         dpar : parameters.ParamReals,
                         rwork : rwork);

               request = (Request) lrequest;
            }


            [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true,
               SetLastError = false)]
            private static extern void dcg_check([In, Out] ref Int32 n,
                                                 [In] IntPtr sol,
                                                 [In] IntPtr rhs,
                                                 [In, Out] ref Int32 request,
                                                 [Out] Int32[] ipar,
                                                 [Out] double[] dpar,
                                                 [Out] IntPtr rwork);

            #endregion

            #region Get

            /// <summary>
            /// Retrieves the number of the current iteration and updates the solution.
            /// </summary>
            /// <param name="parameters">INPUT:The parameters.</param>
            /// <param name="n">INPUT:Sets the size of the problem.</param>
            /// <param name="sol">OUTPUT: Array of size <paramref name="n"/>. Contains the initial approximation to the solution vector. Normally it is equal to 0 or to <paramref name="rhs"/>.</param>
            /// <param name="rhs">OUTPUT: Array of size <paramref name="n"/>. Contains the right-hand side vector.</param>
            /// <param name="rwork">INPUT:Array of size (<paramref name="n"/>,4).</param>
            /// <param name="itercount">OUTPUT: Contains the value of the current iteration number.</param>
            public static  void Get(Parameters parameters,
                                          [In] Int32 n,
                                          [Out] IntPtr sol,
                                          [Out] IntPtr rhs,
                                          [In] IntPtr rwork,
                                          [In, Out] ref Int32 itercount)
            {
               Int32 ln = n;
               Int32 lrequest = 0;
               dcg_get(ref ln,
                       sol,
                       rhs,
                       ref lrequest, // this is not used, see documentation.
                       parameters.ParamIntegers,
                       parameters.ParamReals,
                       rwork,
                       ref itercount);

            }

            [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true,
               SetLastError = false)]
            private static extern  void dcg_get([In, Out] ref Int32 n,
                                                      [Out] IntPtr sol,
                                                      [Out] IntPtr rhs,
                                                      [In, Out] ref Int32 request,
                                                      [In] Int32[] ipar,
                                                      [In] double[] dpar,
                                                      [In] IntPtr rwork,
                                                      [In, Out] ref Int32 itercount);

            #endregion

            #region Run

            /// <summary>
            /// Makes the Conjugate Gradient iterations.
            /// </summary>
            /// <param name="parameters">OUTPUT: The parameters. </param>
            /// <param name="n">INPUT:Sets the size of the problem.</param>
            /// <param name="sol">INPUT: Array of size <paramref name="n"/>. Contains the initial approximation to the solution vector. Normally it is equal to 0 or to b.</param>
            /// <param name="rhs">INPUT: Array of size <paramref name="n"/>. Contains the right-hand side vector.</param>
            /// <param name="rwork">INPUT/OUTPUT:Array of size (n, 4)..</param>
            /// <param name="request">OUTPUT: Gives information about result of the routine.</param>
            public static  void Run(Parameters parameters,
                                          [In] Int32 n,
                                          [In] IntPtr sol,
                                          [In] IntPtr rhs,
                                          [In, Out] IntPtr rwork,
                                          [In, Out] ref Request request)
            {
               Debug.Assert(null != parameters, "The parameter 'parameters' cannot be null.");
               Debug.Assert(null != rwork, "The parameter 'rwork' cannot be null.");
               Debug.Assert(null != sol, "The parameter 'sol' cannot be null.");
               Debug.Assert(null != rhs, "The parameter 'rhs' cannot be null.");

               Int32 ln = n;
               Int32 lrequest = (Int32) request;
               dcg(ref ln,
                   sol,
                   rhs,
                   ref lrequest,
                   parameters.ParamIntegers,
                   parameters.ParamReals,
                   rwork);
               request = (Request) lrequest;
            }


            [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true,
               SetLastError = false)]
            private static extern  void dcg(ref Int32 n,
                                                  IntPtr sol,
                                                  IntPtr rhs,
                                                  ref Int32 request,
                                                  Int32[] ipar,
                                                  double[] dpar,
                                                  IntPtr rwork);

            #endregion

         }

      }

      #endregion

      #region BLAS

      /// <summary>
      /// Routines for BLAS.
      /// </summary>
      public static class BLAS
      {

         #region ?axpy

         /// <summary>
         /// Computes a vector-scalar product and adds the result to a vector, y := a*x + y.
         /// </summary>
         /// <param name="n">Specifies the number of elements in vectors x and y.</param>
         /// <param name="a">Specifies the scalar a.</param>
         /// <param name="x">Array, size at least (1 + (n-1)*abs(<paramref name="incx"/>)).</param>
         /// <param name="incx">Specifies the increment for the elements of x.</param>
         /// <param name="y">Array, size at least (1 + (n-1)*abs(<paramref name="incy"/>)).</param>
         /// <param name="incy">Specifies the increment for the elements of y.</param>
         public static void axpy([In] Int32 n,
                                        [In] double a,
                                        [In] IntPtr x,
                                        [In] Int32 incx,
                                        [Out] IntPtr y,
                                        [In] Int32 incy)
         {
            Debug.Assert(null != x, "The parameter 'x' cannot be null.");
            Debug.Assert(null != y, "The parameter 'y' cannot be null.");

            cblas_daxpy(n,
                        a,
                        x,
                        incx,
                        y,
                        incy);
         }

         [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
         private static extern void cblas_daxpy([In] Int32 n,
                                                       [In] double a,
                                                       [In] IntPtr x,
                                                       [In] Int32 incx,
                                                       [Out] IntPtr y,
                                                       [In] Int32 incy);

         #endregion

         #region nrm2

         /// <summary>
         /// Perform a vector reduction operation defined as res = ||<paramref name="x"/>||,where:
         /// <paramref name="x"/> is a vector,
         /// res is a value containing the Euclidean norm of the elements of <paramref name="x"/>.
         /// </summary>
         /// <param name="n">Specifies the number of elements in vector <paramref name="x"/>.</param>
         /// <param name="x">Array, size at least (1 + (<paramref name="n"/> -1)*abs (<paramref name="incx"/>)).</param>
         /// <param name="incx">Specifies the increment for the elements of <paramref name="x"/>.</param>
         public static  double nrm2([In] Int32 n,
                                          [In] IntPtr x,
                                          [In] Int32 incx)
         {
            return cblas_dnrm2(n,
                               x,
                               incx);
         }

         /// <see cref="nrm2(Int32,IntPtr,Int32)"/>
         public static double nrm2([In] Int32 n,
                                          [In] double[] x,
                                          [In] Int32 incx)
         {
            return nrm2(n,
                        ArrayPointerFactory.Create(x).Handle,
                        incx);
         }

         
         /// <see cref="nrm2(Int32,IntPtr,Int32)"/>
         public static double nrm2([In] int n,
                                   [In] double[] x)
         {
            return nrm2(n, x, 1);
         }
         
         /// <see cref="nrm2(Int32,IntPtr,Int32)"/>
         public  static double nrm2([In] int n,
                                   [In] IntPtr x)
         {
            return nrm2(n, x, 1);
         }

         
         /// <see cref="nrm2(Int32,IntPtr,Int32)"/>
         public static double nrm2([In] double[] x)
         {
            return nrm2(x.Length, x);
         }

      [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
      private static extern double cblas_dnrm2([In] Int32 n,
                                               [In] IntPtr x,
                                               [In] Int32 xld);
         #endregion

         #region gemm

         /// <summary>
         /// Computes a matrix-matrix product with general matrices.
         /// </summary>
         /// <remarks>
         /// The ?gemm routines compute a scalar-matrix-matrix product and add the result to a scalar-matrix product, with general matrices. The operation is defined as
         /// C := alpha*op(A)*op(B) + beta*C,
         /// where:
         /// op(X) is one of op(X) = X, or op(X) = X^T, or op(X) = X^H,
         /// alpha and beta are scalars,
         /// A, B and C are matrices:
         /// op(A) is an <paramref name="m"/>-by-<paramref name="k"/> matrix,
         /// op(B) is a <paramref name="k"/>-by-<paramref name="n"/> matrix,
         /// C is an <paramref name="m"/>-by-<paramref name="n"/> matrix.
         /// </remarks>
         /// <param name="transposeA">Specifies the form of op(A) used in the matrix multiplication:
         /// if <paramref name="transposeA"/> = 'N' or 'n', then op(A) = A;
         /// if <paramref name="transposeA"/> = 'T' or 't', then op(A) = A^T;
         /// </param>
         /// <param name="transposeB">Specifies the form of op(B) used in the matrix multiplication:
         /// if <paramref name="transposeB"/> = 'N' or 'n', then op(B) = B;
         /// if <paramref name="transposeB"/> = 'T' or 't', then op(B) = B^T;
         /// </param>
         /// <param name="m">Specifies the number of rows of the matrix op(A) and of the matrix C. The value of <paramref name="m"/> must be at least zero.</param>
         /// <param name="n"> Specifies the number of columns of the matrix op(B) and the number of columns of the matrix C.</param>
         /// <param name="k">Specifies the number of columns of the matrix op(A) and the number of rows of the matrix op(B).</param>
         /// <param name="alpha">Specifies the scalar <paramref name="alpha"/>.</param>
         /// <param name="a">Array, size <paramref name="lda"/> by ka, where ka is <paramref name="k"/> when transa = 'N' or 'n', and is <paramref name="m"/> otherwise. 
         /// Before entry with transa = 'N' or 'n', the leading <paramref name="m"/>-by-<paramref name="k"/> part of the array a must contain the matrix A, otherwise the 
         /// leading <paramref name="k"/>-by-<paramref name="m"/> part of the array <paramref name="a"/> must contain the matrix A.</param>
         /// <param name="lda">Specifies the leading dimension of <paramref name="a"/> as declared in the calling (sub)program. When transa = 'N' or 'n', 
         /// then <paramref name="lda"/> must be at least max(1, <paramref name="m"/>), otherwise <paramref name="lda"/> must be at least max(1, <paramref name="k"/>).</param>
         /// <param name="b">Array, size <paramref name="ldb"/> by kb, where kb is <paramref name="n"/> when transa = 'N' or 'n', and is <paramref name="k"/> otherwise. 
         /// Before entry with transa = 'N' or 'n', the leading <paramref name="k"/>-by-<paramref name="n"/> part of the array <paramref name="b"/> must contain the matrix B, 
         /// otherwise the leading <paramref name="n"/>-by-<paramref name="k"/> part of the array <paramref name="b"/> must contain the matrix B.</param>
         /// <param name="ldb">When transb = 'N' or 'n', then <paramref name="ldb"/> must be at least max(1, <paramref name="k"/>), otherwise <paramref name="ldb"/> must be at least max(1, <paramref name="n"/>).</param>
         /// <param name="beta">Specifies the scalar <paramref name="beta"/>. When <paramref name="beta"/> is equal to zero, then <paramref name="c"/> need not be set on input.</param>
         /// <param name="c">Array, size <paramref name="ldc"/> by <paramref name="n"/>. Before entry, the leading <paramref name="m"/>-by-<paramref name="n"/> part of the array <paramref name="c"/> must contain the matrix C, 
         /// except when beta is equal to zero, in which case <paramref name="c"/> need not be set on entry.</param>
         /// <param name="ldc">Specifies the leading dimension of <paramref name="c"/> as declared in the calling (sub)program. The value of <paramref name="ldc"/> must be at least max(1, <paramref name="m"/>).</param>
         public static void gemm([In] Int32 transposeA,
                                 [In] Int32 transposeB,
                                 [In] Int32 m,
                                 [In] Int32 n,
                                 [In] Int32 k,
                                 [In] double alpha,
                                 [In] IntPtr a,
                                 [In] Int32 lda,
                                 [In] IntPtr b,
                                 [In] Int32 ldb,
                                 [In] double beta,
                                 [Out] IntPtr c,
                                 [In] Int32 ldc)
         {
            unsafe
            {
               dgemm(&transposeA,
                     &transposeB,
                     &m,
                     &n,
                     &k,
                     &alpha,
                     a,
                     &lda,
                     b,
                     &ldb,
                     &beta,
                     c,
                     &ldc);
            }
         }

         #endregion

         [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
         private static extern unsafe void dgemm([In] Int32* transposeA,
                                                 [In] Int32* transposeB,
                                                 [In] Int32* m,
                                                 [In] Int32* n,
                                                 [In] Int32* p,
                                                 [In] double* a,
                                                 [In] IntPtr x,
                                                 [In] Int32* xld,
                                                 [In] IntPtr y,
                                                 [In] Int32* yld,
                                                 [In] double* s,
                                                 [Out] IntPtr r,
                                                 [In] Int32* rld);







         #region gemm

         /// <summary>
         /// Computes a matrix-matrix product with general matrices.
         /// </summary>
         /// <remarks>
         /// The ?gemm routines compute a scalar-matrix-matrix product and add the result to a scalar-matrix product, with general matrices. The operation is defined as
         /// C := alpha*op(A)*op(B) + beta*C,
         /// where:
         /// op(X) is one of op(X) = X, or op(X) = X^T, or op(X) = X^H,
         /// alpha and beta are scalars,
         /// A, B and C are matrices:
         /// op(A) is an <paramref name="m"/>-by-<paramref name="k"/> matrix,
         /// op(B) is a <paramref name="k"/>-by-<paramref name="n"/> matrix,
         /// C is an <paramref name="m"/>-by-<paramref name="n"/> matrix.
         /// </remarks>
         /// <param name="transposeA">Specifies the form of op(A) used in the matrix multiplication:
         /// if <paramref name="transposeA"/> = 'N' or 'n', then op(A) = A;
         /// if <paramref name="transposeA"/> = 'T' or 't', then op(A) = A^T;
         /// </param>
         /// <param name="transposeB">Specifies the form of op(B) used in the matrix multiplication:
         /// if <paramref name="transposeB"/> = 'N' or 'n', then op(B) = B;
         /// if <paramref name="transposeB"/> = 'T' or 't', then op(B) = B^T;
         /// </param>
         /// <param name="m">Specifies the number of rows of the matrix op(A) and of the matrix C. The value of <paramref name="m"/> must be at least zero.</param>
         /// <param name="n"> Specifies the number of columns of the matrix op(B) and the number of columns of the matrix C.</param>
         /// <param name="k">Specifies the number of columns of the matrix op(A) and the number of rows of the matrix op(B).</param>
         /// <param name="alpha">Specifies the scalar <paramref name="alpha"/>.</param>
         /// <param name="a">Array, size <paramref name="lda"/> by ka, where ka is <paramref name="k"/> when transa = 'N' or 'n', and is <paramref name="m"/> otherwise. 
         /// Before entry with transa = 'N' or 'n', the leading <paramref name="m"/>-by-<paramref name="k"/> part of the array a must contain the matrix A, otherwise the 
         /// leading <paramref name="k"/>-by-<paramref name="m"/> part of the array <paramref name="a"/> must contain the matrix A.</param>
         /// <param name="lda">Specifies the leading dimension of <paramref name="a"/> as declared in the calling (sub)program. When transa = 'N' or 'n', 
         /// then <paramref name="lda"/> must be at least max(1, <paramref name="m"/>), otherwise <paramref name="lda"/> must be at least max(1, <paramref name="k"/>).</param>
         /// <param name="b">Array, size <paramref name="ldb"/> by kb, where kb is <paramref name="n"/> when transa = 'N' or 'n', and is <paramref name="k"/> otherwise. 
         /// Before entry with transa = 'N' or 'n', the leading <paramref name="k"/>-by-<paramref name="n"/> part of the array <paramref name="b"/> must contain the matrix B, 
         /// otherwise the leading <paramref name="n"/>-by-<paramref name="k"/> part of the array <paramref name="b"/> must contain the matrix B.</param>
         /// <param name="ldb">When transb = 'N' or 'n', then <paramref name="ldb"/> must be at least max(1, <paramref name="k"/>), otherwise <paramref name="ldb"/> must be at least max(1, <paramref name="n"/>).</param>
         /// <param name="beta">Specifies the scalar <paramref name="beta"/>. When <paramref name="beta"/> is equal to zero, then <paramref name="c"/> need not be set on input.</param>
         /// <param name="c">Array, size <paramref name="ldc"/> by <paramref name="n"/>. Before entry, the leading <paramref name="m"/>-by-<paramref name="n"/> part of the array <paramref name="c"/> must contain the matrix C, 
         /// except when beta is equal to zero, in which case <paramref name="c"/> need not be set on entry.</param>
         /// <param name="ldc">Specifies the leading dimension of <paramref name="c"/> as declared in the calling (sub)program. The value of <paramref name="ldc"/> must be at least max(1, <paramref name="m"/>).</param>
         public static void gemm([In] Int32 transposeA,
                                 [In] Int32 transposeB,
                                 [In] Int32 m,
                                 [In] Int32 n,
                                 [In] Int32 k,
                                 [In] double alpha,
                                 [In] double[] a,
                                 [In] Int32 lda,
                                 [In] double[]b,
                                 [In] Int32 ldb,
                                 [In] double beta,
                                 [Out] double[] c,
                                 [In] Int32 ldc)
         {
            unsafe
            {
               dgemm(&transposeA,
                     &transposeB,
                     &m,
                     &n,
                     &k,
                     &alpha,
                     a,
                     &lda,
                     b,
                     &ldb,
                     &beta,
                     c,
                     &ldc);
            }
         }

         #endregion

         [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
         private static extern unsafe void dgemm([In] Int32* transposeA,
                                                 [In] Int32* transposeB,
                                                 [In] Int32* m,
                                                 [In] Int32* n,
                                                 [In] Int32* p,
                                                 [In] double* a,
                                                 [In] double[]x,
                                                 [In] Int32* xld,
                                                 [In] double[]y,
                                                 [In] Int32* yld,
                                                 [In] double* s,
                                                 [Out] double[]r,
                                                 [In] Int32* rld);





      }

      #endregion

      #region SPBLAS

      /// <summary>
      /// Routines for the SparseBLAS.
      /// </summary>
      public static class SPBLAS
      {
         #region ?csrgemv

         /// <summary>
         /// General matrix vector product for a sparse matrix with the Compressed Sparse Row format with 0-based indexing.
         /// </summary>
         /// <param name="transa">No transpose ('N' or 'n'), Transpose ('T' or 't')</param>
         /// <param name="m">Number of rows of the matrix A.</param>
         /// <param name="a">Array containing non-zero elements of the matrix A. Its length is equal to the number of non-zero elements in the matrix A.</param>
         /// <param name="ia">Array of length <paramref name="m"/> + 1, containing indices of elements in the array a, such that <paramref name="ia"/>[I] is the index in the array a of the first non-zero element from the row I. The value of the last element <paramref name="ia"/>[m] is equal to the number of non-zeros.</param>
         /// <param name="ja">Array containing the column indices for each non-zero element of the matrix A, its length is equal to the length of the array <paramref name="a"/>. </param>
         /// <param name="x">Input vector.</param>
         /// <param name="y">Output vector.</param>
         [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)
         ]
         public static extern unsafe void mkl_cspblas_dcsrgemv([In] Int32* transa,
                                                               [In] Int32* m,
                                                               [In] IntPtr a,
                                                               [In] IntPtr ia,
                                                               [In] IntPtr ja,
                                                               [In] IntPtr x,
                                                               [Out] IntPtr y);

         /// <summary>
         /// General matrix vector product for a sparse matrix with the Compressed Sparse Row format with 1-based indexing.
         /// </summary>
         /// <param name="transa">No transpose ('N' or 'n'), Transpose ('T' or 't')</param>
         /// <param name="m">Number of rows of the matrix A.</param>
         /// <param name="a">Array containing non-zero elements of the matrix A. Its length is equal to the number of non-zero elements in the matrix A.</param>
         /// <param name="ia">Array of length <paramref name="m"/> + 1, containing indices of elements in the array a, such that <paramref name="ia"/>[I] is the index in the array a of the first non-zero element from the row I. The value of the last element <paramref name="ia"/>[m] is equal to the number of non-zeros.</param>
         /// <param name="ja">Array containing the column indices for each non-zero element of the matrix A, its length is equal to the length of the array <paramref name="a"/>. </param>
         /// <param name="x">Input vector.</param>
         /// <param name="y">Output vector.</param>
         [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
         public static extern unsafe void mkl_dcsrgemv([In] Int32* transa,
                                                       [In] Int32* m,
                                                       [In] IntPtr a,
                                                       [In] IntPtr ia,
                                                       [In] IntPtr ja,
                                                       [In] IntPtr x,
                                                       [Out] IntPtr y);

         #endregion

         #region ?csrtrsv

         /// <summary>
         /// Triangular solvers with simplified interface for a sparse matrix in the CSR format (3-array variation) with one-based indexing.
         /// </summary>
         /// <param name="uplo">Specifies whether the upper or low triangle of the matrix A is used. 
         /// If <paramref name="uplo"/> = 'U' or 'u', then the upper triangle of the matrix A is used. 
         /// If <paramref name="uplo"/> = 'L' or 'l', then the low triangle of the matrix A is used.</param>
         /// <param name="transa">Specifies the system of linear equations.
         /// If <paramref name="transa"/> = 'N' or 'n', then A*y = x.
         /// If <paramref name="transa"/> = 'T' or 't' or 'C' or 'c', then AT*y = x,</param>
         /// <param name="diag">Specifies whether A is unit triangular.
         /// If <paramref name="diag"/> = 'U' or 'u', then A is a unit triangular. 
         /// If <paramref name="diag"/> = 'N' or 'n', then A is not unit triangular.</param>
         /// <param name="m"> Number of rows of the matrix A.</param>
         /// <param name="a">Array containing non-zero elements of the matrix A. Its length is equal to the number of non-zero elements in the matrix A.</param>
         /// <param name="ia">Array of length If <paramref name="m"/> + 1, containing indices of elements in the array <paramref name="a"/>, 
         /// such that <paramref name="ia"/>(i) is the index in the array <paramref name="a"/> of the first non-zero element from the row i. 
         /// The value of the last element <paramref name="ia"/>(<paramref name="m"/> + 1) is equal to the number of non-zeros plus one. </param>
         /// <param name="ja">Array containing the column indices for each non-zero element of the matrix A. Its length is equal to the length of the array a.</param>
         /// <param name="x">Array, size is <paramref name="m"/>. On entry, the array x must contain the vector x.</param>
         /// <param name="y">Array, size at least <paramref name="m"/>. Contains the vector y.</param>
         /// <remarks>
         /// 
         /// The mkl_dcsrtrsv routine solves a system of linear equations with matrix-vector operations for a sparse matrix stored in the CSR format (3 array variation):
         /// A*<paramref name="y"/> = <paramref name="x"/> or A^T * <paramref name="y"/> = <paramref name="x"/>.
         /// where: <paramref name="x"/> and <paramref name="y"/> are vectors, A is a sparse upper or lower triangular matrix with unit or non-unit main diagonal, A^T is the transpose of A.
         /// 
         /// IMPORTANT: This routine supports only one-based indexing of the input arrays.
         /// </remarks>
         [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
         public static extern unsafe void mkl_dcsrtrsv([In] Int32* uplo,
                                                       [In] Int32* transa,
                                                       [In] Int32* diag,
                                                       [In] Int32* m,
                                                       [In] IntPtr a,
                                                       [In] IntPtr ia,
                                                       [In] IntPtr ja,
                                                       [In] IntPtr x,
                                                       [Out] IntPtr y);


         /// <summary>
         /// Triangular solvers with simplified interface for a sparse matrix in the CSR format (3-array variation) with zero-based indexing.
         /// </summary>
         /// <param name="uplo">Specifies whether the upper or low triangle of the matrix A is used. 
         /// If <paramref name="uplo"/> = 'U' or 'u', then the upper triangle of the matrix A is used. 
         /// If <paramref name="uplo"/> = 'L' or 'l', then the low triangle of the matrix A is used.</param>
         /// <param name="transa">Specifies the system of linear equations.
         /// If <paramref name="transa"/> = 'N' or 'n', then A*y = x.
         /// If <paramref name="transa"/> = 'T' or 't' or 'C' or 'c', then AT*y = x,</param>
         /// <param name="diag">Specifies whether A is unit triangular.
         /// If <paramref name="diag"/> = 'U' or 'u', then A is a unit triangular. 
         /// If <paramref name="diag"/> = 'N' or 'n', then A is not unit triangular.</param>
         /// <param name="m"> Number of rows of the matrix A.</param>
         /// <param name="a">Array containing non-zero elements of the matrix A. Its length is equal to the number of non-zero elements in the matrix A.</param>
         /// <param name="ia">Array of length If <paramref name="m"/> + 1, containing indices of elements in the array <paramref name="a"/>, 
         /// such that <paramref name="ia"/>(i) is the index in the array <paramref name="a"/> of the first non-zero element from the row i. 
         /// The value of the last element <paramref name="ia"/>(<paramref name="m"/> + 1) is equal to the number of non-zeros plus one. </param>
         /// <param name="ja">Array containing the column indices for each non-zero element of the matrix A. Its length is equal to the length of the array a.</param>
         /// <param name="x">Array, size is <paramref name="m"/>. On entry, the array x must contain the vector x.</param>
         /// <param name="y">Array, size at least <paramref name="m"/>. Contains the vector y.</param>
         /// <remarks>
         /// 
         /// The mkl_dcsrtrsv routine solves a system of linear equations with matrix-vector operations for a sparse matrix stored in the CSR format (3 array variation):
         /// A*<paramref name="y"/> = <paramref name="x"/> or A^T * <paramref name="y"/> = <paramref name="x"/>.
         /// where: <paramref name="x"/> and <paramref name="y"/> are vectors, A is a sparse upper or lower triangular matrix with unit or non-unit main diagonal, A^T is the transpose of A.
         /// 
         /// IMPORTANT: This routine supports only one-based indexing of the input arrays.
         /// </remarks>
         [DllImport(MklDllName, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)
         ]
         public static extern unsafe void mkl_cspblas_dcsrtrsv([In] Int32* uplo,
                                                               [In] Int32* transa,
                                                               [In] Int32* diag,
                                                               [In] Int32* m,
                                                               [In] IntPtr a,
                                                               [In] IntPtr ia,
                                                               [In] IntPtr ja,
                                                               [In] IntPtr x,
                                                               [Out] IntPtr y);


         #endregion

      }

      #endregion

   }
}
