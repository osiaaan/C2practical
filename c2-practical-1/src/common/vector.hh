#ifndef VECTOR_HH_INCLUDED
#define VECTOR_HH_INCLUDED

#include <iostream>

using namespace std;

/**
     @addtogroup Vector 
     @{
*/


/** \brief 
  class for a variable sized vector of any type (any? why not?)
*/
template <class T>
class Vector 
{
  typedef Vector<T> ThisType;
protected:  
  // constructor
  Vector() : N_(0), data_(0) 
  {
  }

private:
  //! copy constructor
  Vector(const Vector& org) : N_(org.N_) , data_(0)
  {
    // resize vector 
    resize( N_ );
    // copy values 
    for( int i=0; i<N_; ++i ) 
      data_[i] = org.data_[i];
  }

public:
  /** \brief constructor creating vector of length N 
      \param[in] N length of vector 
  */
  Vector(const int N) : data_(0) 
  {
    if( N > 0 )
    {
      resize( N );
    }
    else 
    {
      std::cerr <<"Cannot create Vector of length zero! "<< std::endl;
      assert( false );
      abort();
    }
  }

  /** \brief destructor deleting data */
  ~Vector() 
  {
    delete [] data_;
  }


  /** \brief assignment operator, sets this = v 
      \param[in] other vector that is copied 

      \return reference to this 
  */ 
  Vector& operator=(const Vector& other) 
  {
    // Good coding practice: check for self-assignment in assignment operator
    if (this != &other) 
    {
      // reallocate memory 
      if( size() != other.size() ) 
      {
        resize( other.size() );
      }
      
      for (int i = 0; i < N_; i++)
        data_[i] = other[i];
    }
    return *this;
  }

  /** \brief return size of vector 
      \return size of vector 
  */ 
  size_t size() const 
  {
    assert( N_ >= 0 );
    return N_;
  }

  /** \brief set all entries to zero */
  void clear() 
  {
    assert( data_ );
    for (int i=0; i<N_; ++i) 
    {
      data_[i] = 0;
    }
  }

  /** \brief get value at position i using operator[] (modifiable)
      \param[in] i position of value 
      
      \return const reference to entry i 
  */ 
  const T& operator[](const int i) const {
    checkIndex(i);
    return data_[i];
  }

  /** \brief get value at position i using operator[] (modifiable)
      \param[in] i position of value 
      
      \return reference to entry i 
  */ 
  T& operator[] (const int i) 
  {
    checkIndex(i);
    return data_[i];
  }

  /** \brief scalar product with other vector 
      \param[in] v other vector for evaluating scalar product 

      \return scalar product 
  */ 
  T operator*(const ThisType& v) const 
  {
    if (v.size()!=size()) 
    {
      cerr << "Non compatible dimensions (" << v.size() << " and " << size()
           << ") in Vector::operator*(const Vector& v)" << endl;
      abort();
    }
    T ret = 0; 
    for (int i=0;i<N_; ++i)
    {
      ret += data_[i] * v[i]; 
    }
    return ret;
  }
  
  /** \brief return leak pointer for data (only use for visualization!) 
      \return internal pointer T*   
  */
  T* raw() {
    return data_;
  }

  /** \brief return leak pointer for data (only use for visualization!) 
      \return internal pointer const T*   
  */
  const T* raw() const {
    return data_;
  }

  /** \brief print vector to stream 
      \param s stream that vector is printed to 
  */
  void print(std::ostream& s) const 
  {
    s << "Print Vector: size = "<< size() << " {" <<std::endl; 
    for (int i = 0; i <N_; ++i) 
    {
      s << "v["<< i << "] = " << data_[i] << std::endl; 
    }
    s << "}" << std::endl;
  }
  
  // other useful methods: +=, -=, *=T, addScaled()...
private:
  void checkIndex(int i) const {
    if (i<0 || i>=N_) {
      cerr << "Wrong index " << i << " in Vector::operator[](int i)" 
           << endl;
      abort();
    }
  }

protected:
  /** \brief resize vector with new size, old data is \b lost !!!  
      \param newSize new size of vector 
  */
  void resize( const int newSize ) 
  {
    // set new size 
    N_ = newSize;
    
    // reallocate memory 
    delete [] data_;
    data_ = new T[N_];
    
    // init with zero 
    clear();
  }
  
  //! size of vector
  int N_;       
  //! data array
  T* data_;    
};

/** \brief overloaded operator << for printing vector 
 *
    \param s stream that the vector is printed to 
    \param v vector to print 
*/
template <class T>
std::ostream& operator << (std::ostream& s, const Vector<T>& v) 
{
  v.print(s);
  return s;
}

/** 
   @}
*/

#endif
