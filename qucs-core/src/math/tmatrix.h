/*
 * tmatrix.h - simple matrix template class definitions
 *
 * Copyright (C) 2004, 2005, 2006 Stefan Jahn <stefan@lkcc.org>
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this package; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.  
 *
 * $Id$
 *
 */

#ifndef __TMATRIX_H__
#define __TMATRIX_H__

#include <Eigen/Core>
#include <Eigen/Dense>

#include <assert.h>



template <typename nr_type_t> class tvector;
template <typename nr_type_t> class tmatrix;
template <typename T>  tmatrix<T> operator * (const tmatrix<T> &a, const tmatrix<T> & b);
template <typename T>  tmatrix<T> operator * (const T &a, const tmatrix<T> & b);
template <typename nr_type_t> tvector<nr_type_t> operator * (const tmatrix<nr_type_t> &a, const tvector<nr_type_t> &b);

/*! \brief Compatibility class arround Eigen 
    \todo: die
*/
template <typename nr_type_t>
class tmatrix
{
 public:
  /*! \brief default constructor */
  tmatrix () : m() {};

  /*! \brief Create a square matrix */
  // tmatrix (const int i): m(i,i) {};

  /* \brief General constructor case */
  tmatrix (const int line, const int col) : m(line,col) {};

  /* default copy constructor */
  /* tmatrix (const tmatrix &); */

  /* \brief Creates a tmatrix from an eigen matrix */
  template <typename Derived> tmatrix(const Eigen::EigenBase<Derived>& a) {
    this->m = a;
  };
  
  /* \brief conversion operator to eigen class */
  template <typename Derived> operator Eigen::EigenBase<Derived>() {
    return this->m;
  };

  /*\brief assignment copy constructor 
    \param[in] a: matrix to assign
  */
  const tmatrix& operator = (const tmatrix & a) {
    if(&a != this) {
      this->m = a.m;
    }
    return *this;
  }

  /*! easy accessor operators (read only) */
  nr_type_t  operator () (const int r, const int c) const {
    return m(r,c);
  }
  
  /*! easy accessor operators  */
  nr_type_t& operator () (const int r, const int c) {
    return m(r,c);
  }

  /*! return number of columns */
  int cols() const {
    return this->m.cols(); 
  }

  /*! return number of rows */
  int rows() const {
    return this->m.rows(); 
  }

  /* default destructor thanks eigen */
  /* ~tmatrix (); */
  
  /*!\brief fill matrix */
  void setConstant (const nr_type_t &v) {
     this->m.setConstant(v);
  }
  
  /*! \brief return raw array */
  nr_type_t * const data (void) const { 
    return this->m.data(); 
  }
  
  /*!\brief The function returns the given row in a tvector 
    \todo simplifie when tvector is converted
  */
  tvector<nr_type_t> getRow (const int r) const {
    const int cols = this.cols();
    const int rows = this.rows();
    
    assert (r >= 0 && r < rows);
    
    tvector<nr_type_t> res (cols);
    for(int c = 0; c < cols; c++)
      res(c) = this->m(r,c);
   
    return res;
  }

  /*! Puts the given tvector into the given row of the tmatrix instance */
  void setRow (const int r, const tvector<nr_type_t> & v) {
    const int cols = this.cols();
    const int rows = this.rows();   

    assert (r >= 0 && r < rows && v.getSize () == cols);

    for(int c = 0; c < cols; c++)
      this->m(r,c) = v(c);
  }
  
  /*! The function returns the given column in a tvector */
  tvector<nr_type_t> getCol (const int c) const {
    const int cols = this.cols();
    const int rows = this.rows();
    assert (c >= 0 && c < cols);

    tvector<nr_type_t> res (rows);
    for(int r = 0; r < rows; r++)
      res(r) = this->m(r,c);
    return res;
  }

  /*! Puts the given tvector into the given column of the tmatrix instance. */
  void setCol (const int c, const tvector<nr_type_t> &v) {
    const int cols = this.cols();
    const int rows = this.rows();
    assert (c >= 0 && c < cols);

    for(int r = 0; r < rows; r++)
      this->m(r,c) = v(c);
  }

  /*!\brief The function swaps the given rows with each other.
    \param[in] r1 source row
    \param[in] r2 destination row
  */
  void exchangeRows (const int r1, const int r2) {
    this->m.row(r1).swap(this->m.row(r2));
  }

  /*!\brief The function swaps the given column with each other.
    \param[in] c1 source column
    \param[in] c2 destination column
  */
  void exchangeCols (const int c1, const int c2) {
    this->m.col(c1).swap(this->m.col(c2));
  }
  
  /*! \transpose in place a matrix */
  void transposeInPlace (void) {
    this->m.transposeInPlace();
  }
  
  /*! Checks validity of matrix*/
  bool  isFinite (void) const {
    const int rows = this->rows();
    const int cols = this->cols();
    for (int i = 0; i < rows; i++)
      for(int j=0;j < cols; j++)
	if (!finite (real ((*this)(i,j)))) 
	  return false;
    return true;
  };
  
  void print (bool realonly = false);


 /*!\brief return L inf norm of norm-1 vector of rowise */
 nr_double_t infnorm () {
   return ((this->m.rowwise()).template lpNorm<1>()).template lpNorm<Eigen::Infinity>();
 }

 /* \brief return condition number 
    \note slow For debugging purposes only
 */
 nr_double_t condition () {
   return this->infnorm () * inverse(this).infnorm();
 }

  template<typename T> friend tmatrix<T> operator * (const tmatrix<T> &a, const tmatrix<T> & b);
  template<typename T> friend tmatrix<T> operator * (const T &a, const tmatrix<T> & b);
  template<typename T> friend tvector<T> operator * (const tmatrix<T> &a, const tvector<T> &b);
  template<typename T> friend tvector<T> operator * (const tvector<T> &a, const tmatrix<T> &b);
  template<typename T> friend tmatrix<T> inverse (const tmatrix<T> &a);
    
  static tmatrix Identity(const int n, const int m) {
    tmatrix<nr_type_t> ret;
    ret.m = Eigen::Matrix<nr_type_t, Eigen::Dynamic, Eigen::Dynamic >::Identity(n,n);
    return ret;
  }
 
  // intrinsic operators
  tmatrix<nr_type_t> operator += (const tmatrix<nr_type_t> & t) {
    this->m+=t.m;
    return *this;
  }

  tmatrix<nr_type_t> operator -= (const tmatrix<nr_type_t> &t) {
    this->m-= t.m;
    return *this;
  }

  // intrinsic operators
  tmatrix<nr_type_t> operator *= (const nr_type_t & t) {
    this->m*=t;
    return *this;
  }

 private:
  /*! Matrix data */
  Eigen::Matrix<nr_type_t,Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> m;
};

/*! Multiplication of two matrix */
template <typename nr_type_t>
tmatrix<nr_type_t> operator * (const tmatrix<nr_type_t> &a, const tmatrix<nr_type_t>& b) {
  assert (a.cols () == b.rows ());

  tmatrix<nr_type_t> ret;
  ret.m = a.m * b.m;
  
  return ret;
}

/*! Multiplication of scalar and matrix */
template <typename nr_type_t>
tmatrix<nr_type_t> operator * (const nr_type_t &a, const tmatrix<nr_type_t>& b) {
  return a * b.m;
}


/*! Multiplication of matrix and vector
    \todo should die 
*/
template <class nr_type_t>
tvector<nr_type_t> operator * (const tmatrix<nr_type_t> &a, const tvector<nr_type_t> &b) {
  assert (a.cols () == b.size ());
  int r, c, n = a.cols ();
  nr_type_t z;
  tvector<nr_type_t> res (n);

  for (r = 0; r < n; r++) {
    for (c = 0, z = 0; c < n; c++) z += a(r, c) * b(c);
    res.set (r, z);
  }
  return res;
}

template <class nr_type_t>
tvector<nr_type_t> operator * (const tvector<nr_type_t> &a, const tmatrix<nr_type_t> &b) {
  assert (a.size () == b.rows ());
  int r, c, n = b.rows ();
  nr_type_t z;
  tvector<nr_type_t> res (n);

  for (c = 0; c < n; c++) {
    for (r = 0, z = 0; r < n; r++) z += a(r) * b(r, c);
    res.set (c, z);
  }
  return res;
}


#include "tmatrix.cpp"

#endif /* __TMATRIX_H__ */
