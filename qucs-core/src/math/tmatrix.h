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

#include "tvector.h"

#include <assert.h>

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
 private:
  /*! Matrix data */
  Eigen::Matrix<nr_type_t,Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> m;
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

  auto array(void) -> decltype(this->m.array()) {
    return m.array();
  }

  auto row(const unsigned int i) -> decltype(this->m.row(i)) {
    return this->m.row(i);
  }

  auto col(const unsigned int i) -> decltype(this->m.col(i)) {
    return this->m.col(i);
  }
  
  auto rowwise() -> decltype(this->m.rowwise()) {
    return m.rowwise();
  }

  template<typename T> friend tmatrix<T> operator * (const tmatrix<T> &a, const tmatrix<T> & b);
  template<typename T> friend tmatrix<T> operator * (const T &a, const tmatrix<T> & b);
  template<typename T> friend tvector<T> operator * (const tmatrix<T> &a, const tvector<T> &b);
  template<typename T> friend tvector<T> operator * (const tvector<T> &a, const tmatrix<T> &b);

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

  // Compute inverse matrix of the given matrix by Gauss-Jordan elimination.
  tmatrix<nr_type_t> inverse () {
    auto & a = *this;
    nr_double_t MaxPivot;
    nr_type_t f;
    tmatrix<nr_type_t> b;
    tmatrix<nr_type_t> e;
    int i, c, r, pivot, n = a.cols ();
  
    // create temporary matrix and the result matrix
    b = tmatrix<nr_type_t> (a);
    e = tmatrix<nr_type_t>().Identity(n,n);
    
    // create the eye matrix in 'b' and the result in 'e'
    for (i = 0; i < n; i++) {
      // find maximum column value for pivoting
      for (MaxPivot = 0, pivot = r = i; r < n; r++) {
	if (abs (b(r, i)) > MaxPivot) {
	  MaxPivot = abs (b(r, i));
	  pivot = r;
	}
      }
      // exchange rows if necessary
      assert (MaxPivot != 0); // singular matrix
      if (i != pivot) {
	b.row(i).swap(b.row(pivot));
	e.row(i).swap(e.row(pivot));
      }
      
      // compute current row
      f = b (i, i);
      for (c = 0; c < n; c++) {
	b(i, c) =  b (i, c) / f;
	e(i, c) =  e (i, c) / f;
      }
      
      // compute new rows and columns
      for (r = 0; r < n; r++) {
	if (r != i) {
	  f = b (r, i);
	  for (c = 0; c < n; c++) {
	    b(r, c) =  b (r, c) - f * b (i, c);
	    e(r, c) =  e (r, c) - f * e (i, c);
	  }
	}
      }
    }
    return e;
  }
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
  tvector<nr_type_t> res = tvector<nr_type_t>::Zero (n,1);

  for (r = 0; r < n; r++) {
    for (c = 0, z = 0; c < n; c++) z += a(r, c) * b(c);
    res(r) = z;
  }
  return res;
}

template <class nr_type_t>
tvector<nr_type_t> operator * (const tvector<nr_type_t> &a, const tmatrix<nr_type_t> &b) {
  assert (a.size () == b.rows ());
  int r, c, n = b.rows ();
  nr_type_t z;
  tvector<nr_type_t> res = tvector<nr_type_t>::Zero (n,1);

  for (c = 0; c < n; c++) {
    for (r = 0, z = 0; r < n; r++) z += a(r) * b(r, c);
    res(c) = z;
  }
  return res;
}


#include "tmatrix.cpp"

#endif /* __TMATRIX_H__ */
