#ifndef __header_MatrixHandle_hpp__
#define __header_MatrixHandle_hpp__

template <typename T> class MatrixHandleReadOnly
{
 protected:

  I 		m_numRows;
  I 		m_numCols;
  I 		m_ld;
  const T * 	m_roPtr;
  
 public:
  
  MatrixHandleReadOnly(const I 	numRows_,
		       const I 	numCols_,
		       const T*	roPtr_,
		       const I 	ld_)
  {
    this->m_numRows	= numRows_;
    this->m_numCols	= numCols_;
    this->m_ld 		= ld_;
    this->m_roPtr 	= roPtr_;
  };

  I GetNumRows() const
  {
    return this->m_numRows;
  };

  I GetLd() const
  {
    return this->m_ld;
  };
  
  I GetNumColumns() 	const
  {
    return this->m_numCols;
  };
  
  T Get(const I i_,const I j_) const
  {
    return this->m_roPtr[this->m_ld*j_+i_];
  };

  virtual ~MatrixHandleReadOnly()
    {
      this->m_numRows 	= 0;
      this->m_numCols 	= 0;
      this->m_ld 	= 0;
      this->m_roPtr 	= NULL;
    };

};

template <typename T> class MatrixHandle : public MatrixHandleReadOnly<T>
{
 protected:

  T * 	m_ptr;
  
 public:
  
  MatrixHandle(const I 	numRows_,
	       const I 	numCols_,
	       T*	ptr_,
	       const I 	ld_) : MatrixHandleReadOnly<T>(numRows_,
						       numCols_,
						       ptr_,
						       ld_)
  {
    this->m_ptr	= ptr_;
  };

  virtual ~MatrixHandle()
    {
      this->m_ptr 	= NULL;
    };


  void Set(const I i_,
	   const I j_,
	   const T x_) 
  {
    this->m_ptr[this->m_ld*j_+i_] = x_;
  };

  
};


template <class T> std::ostream& operator<<(std::ostream&s, const MatrixHandleReadOnly<T> & b)
{
  const I n = b.GetNumRows();
  const I m = b.GetNumColumns();
  s << n << " " << m << endl;
  for (I i=0;i<n;++i,s << endl)
    for (I j=0;j<m;++j)
      s << " " << b.Get(i,j);
  return s;
}


#endif
