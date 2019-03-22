#ifndef SPARSEMATRIX_HH_INCLUDED 
#define SPARSEMATRIX_HH_INCLUDED 

#include <cstdlib>
#include <iostream>

#include "../common/cghs.hh"

/** 
    @addtogroup Matrix 
    @{
*/

/** \brief matrix implementation that stores entries in 
    compressed-row storage (CRS) 
*/ 
template <class T>
class SparseMatrix {
 public:

  /** \brief constructor 
      \param rows number of rows 
      \param nonZero maximal number of non-zeros per row 
  */
  SparseMatrix(const int rows, const int nonZero) 
    : rows_(rows), nonZero_(nonZero) 
  {
    cols_=new int[rows_*nonZero_];
    vals_=new T[rows_*nonZero_];
    
    for (int i=0;i<rows_*nonZero_;i++) 
    {
     cols_[i] = -1;
     vals_[i] = 0.;
    }
  }

  /** \brief destructor */
  ~SparseMatrix() 
  {
    delete [] vals_;
    delete [] cols_;
  }

  /** \brief add value to matrix entry 
      \param row row of matrix entry 
      \param column column of matrix entry 
      \param value value to add 
  */
  void add(const int row , const int column, const T &value) 
  {
    assert(row < rows_);
    assert(value==value);
    const int k = colIndex(row, column);
    vals_[row * nonZero_ + k] += value;
  }
  
  /** \brief set matrix entry to certain value 
      \param row row of matrix entry 
      \param column column of matrix entry 
      \param value value to set to matrix entry 
  */
  void set(const int row, const int column ,const T value) 
  {
    assert(row < rows_);
    assert(value==value);
    const int k=colIndex(row, column);
    vals_[row * nonZero_ + k] = value;
  }
  
  /** \brief set diagonal entry to certain value and remove all other
      column entries for this row 
      
      \param row row of matrix entry 
      \param value value to set to matrix entry (default value is 1)
  */
  void setDiag(const int row, const T value = 1 ) 
  {
    assert(row < rows_);
    assert(value==value);
    cols_[row * nonZero_ ] = row;
    vals_[row * nonZero_ ] = value;

    for (int k=1; k < nonZero_; ++k) 
    {
      cols_[row * nonZero_ + k] = -1;
      vals_[row * nonZero_ + k] = 0;
    }
  }

  /** \brief return number of rows 
      \return number of rows 
  */
  int row() const 
  {
    return rows_;
  }
  
  /** \brief print matrix to standard out stream, for debugging only */
  void print() const 
  {
    for (int r=0; r<rows_; ++r) 
    {
      int k=0;
      int c = getcol(r,k);
      std::cout << "Row " << r << " = "; 
      while (k<maxcol() && c!=-1) 
      {
        std::cout << "(" << c << " : " << getval(r,k) << ") ";
        ++k;
        c=getcol(r,k);
      }
      std::cout << std::endl;
    }
  }
  
  /** \brief matrix-vector multiplication, i.e. \f$ w = A \, v \f$ 
      \param[in]  v vector \f$v\f$ to multiply 
      \param[out] w result vector \f$w\f$  
   */
  void mult(const T* v, T* w ) const 
  {
    for (int r=0; r<row(); ++r) 
    {
      // initialize vector value  
      w[r] = 0;

      for(int k=0; k<nonZero_; ++k) 
      {
        int c = getcol(r,k);
        if( c == -1 ) continue ;
          assert(c>=0);
        assert( v[c] == v[c] );
        // multiply 
        w[r] += getval(r,k) * v[c];
        assert( w[r] == w[r] );
      }
    }
  }
 private:
  int getcol(const int row, const int k) const 
  {
    return cols_[ row * nonZero_ + k];
  }
  
  T getval(const int row, const int k) const 
  {
    assert( vals_[row*nonZero_+k] == vals_[row*nonZero_+k]);
    return vals_[row*nonZero_+k];
  }
  
  int maxcol() const  
  {
    return nonZero_;
  }

  // index of column 
  int colIndex(const int r, const int c) 
  {
    int ret=0;
    while (cols_[ r* nonZero_ + ret] != -1 &&
           cols_[r*nonZero_+ret]!=c) 
    {
      ++ret;
      if (ret == nonZero_) 
      {
       std::cerr << "SparseMatrix: zu wenig Eintraege pro Zeile!" << std::endl;
       abort();
      }
    }
    if (cols_[r*nonZero_+ret] == -1)
      cols_[ r * nonZero_ + ret] = c;
    return ret;
  }
  
  // number of rows 
  int rows_;
  // maximal number of non-zeros per column 
  int nonZero_;
  // pointer for column vector 
  int* cols_;
  // pointer to data vector 
  T* vals_;
};

/** 
  @}
*/

#endif
