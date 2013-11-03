/*
 * tvector.cpp - simple vector template class implementation
 *
 * Copyright (C) 2004, 2005, 2006, 2008 Stefan Jahn <stefan@lkcc.org>
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

#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "compat.h"
#include "complex.h"
#include "tvector.h"

// The function swaps the given rows with each other.
template <class nr_type_t>
void tvector<nr_type_t>::exchangeRows (int r1, int r2) {
  assert (r1 >= 0 && r2 >= 0 && r1 < this->size() && r2 < this->size());
  nr_type_t s = (*this)(r1);
  (*this)(r1) = (*this)(r2);
  (*this)(r2) = s;
}

// Addition.
template <class nr_type_t>
tvector<nr_type_t> operator + (tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert (a.size () == b.size ());
  int n = a.size ();
  tvector<nr_type_t> res = tvector<nr_type_t>::Zero(n,1);
  for (int i = 0; i < n; i++) 
    res(i) = a(i) + b(i);
  return res;
}

// Intrinsic vector addition.
template <class nr_type_t>
tvector<nr_type_t> tvector<nr_type_t>::operator += (tvector<nr_type_t> a) {
  assert (a.size () == this->size());
  for (int i = 0; i < this->size(); i++)
    (*this)(i) += a(i);
  return *this;
}

// Subtraction.
template <class nr_type_t>
tvector<nr_type_t> operator - (tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert (a.size () == b.size ());
  int n = a.size ();
  tvector<nr_type_t> res;
  res = tvector<nr_type_t>::Zero (n,1);
  for (int i = 0; i < n; i++) 
    res(i) =  a(i) - b(i);
  return res;
}

// Intrinsic vector substration.
template <class nr_type_t>
tvector<nr_type_t> tvector<nr_type_t>::operator -= (tvector<nr_type_t> a) {
  assert (a.size () == this->size());
  for (int i = 0; i < this->size(); i++)
    (*this)(i) += a(i);
  return *this;
}

// Intrinsic scalar multiplication.
template <class nr_type_t>
tvector<nr_type_t> tvector<nr_type_t>::operator *= (nr_double_t s) {
  for (int i = 0; i < this->size(); i++) 
    (*this)(i) *= s;
  return *this;
}

// Scalar multiplication.
template <class nr_type_t>
tvector<nr_type_t> operator * (const nr_type_t &s, const tvector<nr_type_t> &a) {
  return tvector<nr_type_t>(s*a.v);
}

template <class nr_type_t>
tvector<nr_type_t> operator * (const tvector<nr_type_t> a, const nr_double_t s) {
  return s * a;
}

// Vector negation.
template <class nr_type_t>
tvector<nr_type_t> operator - (tvector<nr_type_t> a) {
  int n = a.size ();
  tvector<nr_type_t> res = tvector<nr_type_t>::Zero(n,1);
  for (int i = 0; i < n; i++) 
    res(i) =-a(i);
  return res;
}

// Vector less comparison.
template <class nr_type_t>
bool operator < (tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert (a.size () == b.size ());
  int n = a.size ();
  for (int i = 0; i < n; i++) if (a(i) >= b(i)) return false;
  return true;
}

// Vector greater comparison.
template <class nr_type_t>
bool operator > (tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert (a.size () == b.size ());
  int n = a.size ();
  for (int i = 0; i < n; i++) if (a(i) <= b(i)) return false;
  return true;
}

// Scalar addition.
template <class nr_type_t>
tvector<nr_type_t> operator + (nr_type_t s, tvector<nr_type_t> a) {
  int n = a.size ();
  tvector<nr_type_t> res (n);
  for (int i = 0; i < n; i++) 
    res(i) = s + a(i);
  return res;
}

template <class nr_type_t>
tvector<nr_type_t> operator + (tvector<nr_type_t> a, nr_type_t s) {
  return s + a;
}

// Mean square norm.
template <class nr_type_t>
nr_double_t norm (tvector<nr_type_t> a) {
#if 0
  nr_double_t k = 0;
  for (int i = 0; i < a.size (); i++) k += norm (a(i));
  return n;
#else
  nr_double_t scale = 0, n = 1, x, ax;
  for (int i = 0; i < a.size (); i++) {
    if ((x = real (a (i))) != 0) {
      ax = fabs (x);
      if (scale < ax) {
	x = scale / ax;
	n = 1 + n * x * x;
	scale = ax;
      }
      else {
	x = ax / scale;
	n += x * x;
      }
    }
    if ((x = imag (a (i))) != 0) {
      ax = fabs (x);
      if (scale < ax) {
	x = scale / ax;
	n = 1 + n * x * x;
	scale = ax;
      }
      else {
	x = ax / scale;
	n += x * x;
      }
    }
  }
  return scale * scale * n;
#endif
}

// Maximum norm.
template <class nr_type_t>
nr_double_t maxnorm (tvector<nr_type_t> a) {
  nr_double_t nMax = 0, n;
  for (int i = 0; i < a.size (); i++) {
    n = norm (a(i));
    if (n > nMax) nMax = n;
  }
  return nMax;
}

// The functions reorders the vector according to the given index array.
template <class nr_type_t>
void tvector<nr_type_t>::reorder (int * idx) {
  tvector<nr_type_t> old = *this;
  for (int i = 0; i < this->size(); i++) 
    (*this)(i) = old(idx[i]);
}

#ifdef DEBUG
// Debug function: Prints the vector object.
template <class nr_type_t>
void tvector<nr_type_t>::print (bool realonly) {
  for (int r = 0; r < size; r++) {
    if (realonly)
      fprintf (stderr, "%+.2e\n", (double) real (this(r)));
    else
      fprintf (stderr, "%+.2e%+.2ei\n", (double) real (this(r)),
	       (double) imag (this(r)));
  }
  fprintf (stderr, ";;;\n"); 
}
#endif /* DEBUG */
