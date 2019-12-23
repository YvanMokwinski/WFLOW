#ifndef __LinearSolver_DiagonalMatrix_hpp__
#define __LinearSolver_DiagonalMatrix_hpp__

namespace LinearSolver
{

   /// <summary>
   /// Implementation of a diagonal matrix.
   /// </summary>
   class DiagonalMatrix
   {
   protected :

     /// <summary>
     /// The size of the diagonal matrix.
     /// </summary>
     I m_size;

     /// <summary>
     /// The values of the diagonal matrix.
     /// </summary>
     pR m_values;
     
   public:

     /// <summary>
     /// Constructor.
     /// </summary>
     /// <param name="size">The size of the diagonal matrix.</param>
     DiagonalMatrix(const I size)
      {
	this->m_size = size;
	this->m_values = (pR)malloc(sizeof(R)*size);
      };
     
     /// <summary>
     /// Constructor.
     /// </summary>
     /// <param name="size">The size of the diagonal matrix.</param>
     virtual ~DiagonalMatrix()
     {
       if (NULL != this->m_values)
	 {
	   free(this->m_values);
	   this->m_values = NULL;
	 }
       this->m_size = 0;
     };
     
     
     /// <summary>
     /// Get/Set the value of the diagonal coefficient.
     /// </summary>
     /// <param name="index">The index of the diagonal coefficient.</param>
     /// <returns>The value of the diagonal coefficient.</returns>
     R operator[](const I index) const 
     {
       return m_values[index];
     };

     /// <summary>
     /// Compute the inverse.
     /// </summary>
     /// <param name="minAbsoluteValue">Threshold to apply on the absolute value of the diagonal coefficient,
     /// the default value is <see cref="MathAndStats.MachinePrecisionDouble"/>.</param>
     void Inverse(const R minAbsoluteValue = ((R)1.1e-16))
     {
       for (I index = 0; index < this->m_size; ++index)
         {
	   R v = this->m_values[index];
	   
	   if (((v < ((R)0.0)) ? -v : v) < minAbsoluteValue)
	     {
               v = (v < ((R)0.0)) ? -minAbsoluteValue : minAbsoluteValue;
	     }
	   
	   this->m_values[index] = ((R)1.0) / v;
         }
     };
     
     /// <summary>
     /// Perform the matrix-vector product.
     /// </summary>
     /// <param name="transpose">No transpose ('N' or 'n'), Transpose ('T' or 't')</param>
     /// <param name="y">The output vector.</param>
     /// <param name="x">The input vector.</param>
     void Amux(const char* transpose,
	       pR 	y,
	       cst_pR 	x)
     {
       if ((transpose[0] == 'N') || (transpose[0] == 'n') || (transpose[0] == 'T') || (transpose[0] == 't'))
         {
	   for (I index = 0; index < this->m_size; ++index)
	     {
               y[index] = x[index] * this->m_values[index];
	     }
         }       
     };
     

   };

};

#endif
