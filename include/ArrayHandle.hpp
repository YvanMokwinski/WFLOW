#ifndef __header_ArrayHandle_hpp__
#define __header_ArrayHandle_hpp__


template <typename T> class ArrayHandleReadOnly
{
 protected:
  
  I 		m_size;
  const T * 	m_roPtr;
  I 		m_ld;

 public:

  I GetSize() const
  {
    return this->m_size;
  };

  I GetLd() const
  {
    return this->m_ld;
  };

  T Get(const I i) const
  {
    return this->m_roPtr[this->m_ld * i];
  };

  const T * GetPtr() 
  {
    return this->m_roPtr;
  };

  
 public:
  
  ArrayHandleReadOnly(const  	I 	size_,
		      const 	T * 	roPtr_,
		      const 	I 	ld_ = 1)
  {
    this->m_size 	= size_;
    this->m_roPtr 	= roPtr_;
    this->m_ld 		= ld_;
  };
  
  virtual ~ArrayHandleReadOnly()
  {
    this->m_size 	= 0;
    this->m_roPtr  	= NULL;
    this->m_ld 		= 0;
  };
  
};

template <typename T> class ArrayHandle : public ArrayHandleReadOnly<T>
{

 protected:
  T * m_ptr;
  
 public:

  virtual ~ArrayHandle()
  {
    this->m_ptr = NULL;
  };
  
  ArrayHandle(const  I 	size_,
	      T * 	ptr_,
	      const I 	ld_ = 1) : ArrayHandleReadOnly<T>(size_,
							  ptr_,
							  ld_)
  {
    this->m_ptr = ptr_;
  };
  
};



template <class T> std::ostream& operator<<(std::ostream&s, const ArrayHandleReadOnly<T> & b)
{
  const I size = b.GetSize();
  s << size;
  for (I i=0;i<size;++i)
    s << " " << b.Get(i);
  return s;
}



#endif
