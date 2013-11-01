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

// Constructor creates an unnamed instance of the tvector class.
template <class nr_type_t>
tvector<nr_type_t>::tvector () {
  size_ = 0;
  data = NULL;
}

/* Constructor creates an unnamed instance of the tvector class with a
   certain length. */
template <class nr_type_t>
tvector<nr_type_t>::tvector (int s)  {
  size_ = s;
  if (s > 0) {
    data = new nr_type_t[s];
    memset (data, 0, sizeof (nr_type_t) * s);
  }
  else data = NULL;
}

/* The copy constructor creates a new instance based on the given
   tvector object. */
template <class nr_type_t>
tvector<nr_type_t>::tvector (const tvector & v) {
  size_ = v.size_;
  data = NULL;

  // copy tvector elements
  if (this->size_ > 0) {
    data = new nr_type_t[size_];
    memcpy (data, v.data, sizeof (nr_type_t) * size_);
  }
}

/* The assignment copy constructor creates a new instance based on the
   given tvector object. */
template <class nr_type_t>
const tvector<nr_type_t>&
tvector<nr_type_t>::operator=(const tvector<nr_type_t> & v) {
  if (&v != this) {
    size_ = v.size_;
    if (data) { delete[] data; data = NULL; }
    if (size_ > 0) {
      data = new nr_type_t[size_];
      memcpy (data, v.data, sizeof (nr_type_t) * size_);
    }
  }
  return *this;
}

// Destructor deletes a tvector object.
template <class nr_type_t>
tvector<nr_type_t>::~tvector () {
  if (data) delete[] data;
}

// Sets all the tvector elements to the given value.
template <class nr_type_t>
void tvector<nr_type_t>::setConstant (nr_type_t z) {
  for (int i = 0; i < size_; i++) data[i] = z;
}

// Sets the specified tvector elements to the given value.
template <class nr_type_t>
void tvector<nr_type_t>::set (nr_type_t z, int start, int stop) {
  for (int i = start; i < stop; i++) data[i] = z;
}


// Copies the specified elements from the given tvector.
template <class nr_type_t>
void tvector<nr_type_t>::set (tvector<nr_type_t> a, int start, int stop) {
  for (int i = start; i < stop; i++) data[i] = a(i);
} 

// The function swaps the given rows with each other.
template <class nr_type_t>
void tvector<nr_type_t>::exchangeRows (int r1, int r2) {
  assert (r1 >= 0 && r2 >= 0 && r1 < this->size() && r2 < this->size());
  nr_type_t s = data[r1];
  data[r1] = data[r2];
  data[r2] = s;
}

// Addition.
template <class nr_type_t>
tvector<nr_type_t> operator + (tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert (a.size () == b.size ());
  int n = a.size ();
  tvector<nr_type_t> res (n);
  for (int i = 0; i < n; i++) 
    res(i) = a(i) + b(i);
  return res;
}

// Intrinsic vector addition.
template <class nr_type_t>
tvector<nr_type_t> tvector<nr_type_t>::operator += (tvector<nr_type_t> a) {
  assert (a.size () == this->size());
  nr_type_t * src = a.getData ();
  nr_type_t * dst = data;
  for (int i = 0; i < this->size(); i++) *dst++ += *src++;
  return *this;
}

// Subtraction.
template <class nr_type_t>
tvector<nr_type_t> operator - (tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert (a.size () == b.size ());
  int n = a.size ();
  tvector<nr_type_t> res (n);
  for (int i = 0; i < n; i++) 
    res(i) =  a(i) - b(i);
  return res;
}

// Intrinsic vector substration.
template <class nr_type_t>
tvector<nr_type_t> tvector<nr_type_t>::operator -= (tvector<nr_type_t> a) {
  assert (a.size () == this->size());
  nr_type_t * src = a.getData ();
  nr_type_t * dst = data;
  for (int i = 0; i < this->size(); i++) *dst++ -= *src++;
  return *this;
}

// Intrinsic scalar multiplication.
template <class nr_type_t>
tvector<nr_type_t> tvector<nr_type_t>::operator *= (nr_double_t s) {
  nr_type_t * dst = data;
  for (int i = 0; i < this->size(); i++) *dst++ *= s;
  return *this;
}

// Intrinsic scalar division.
template <class nr_type_t>
tvector<nr_type_t> tvector<nr_type_t>::operator /= (nr_double_t s) {
  nr_type_t * dst = data;
  for (int i = 0; i < size; i++) *dst++ /= s;
  return *this;
}

// Scalar multiplication.
template <class nr_type_t>
tvector<nr_type_t> operator * (nr_double_t s, tvector<nr_type_t> a) {
  int n = a.size ();
  tvector<nr_type_t> res (n);
  for (int i = 0; i < n; i++) 
    res(i) = s * a(i);
  return res;
}

template <class nr_type_t>
tvector<nr_type_t> operator * (tvector<nr_type_t> a, nr_double_t s) {
  return s * a;
}

// Vector multiplication (element by element).
template <class nr_type_t>
tvector<nr_type_t> operator * (tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert (a.size () == b.size ());
  int n = a.size ();
  tvector<nr_type_t> res (n);
  for (int i = 0; i < n; i++) 
    res(i) = a(i) * b(i);
  return res;
}

// Computes the scalar product of two vectors.
template <class nr_type_t>
nr_type_t scalar (tvector<nr_type_t> a, tvector<nr_type_t> b) {
  assert (a.size () == b.size ());
  nr_type_t n = 0;
  for (int i = 0; i < a.size (); i++) 
    n += a(i) * b(i);
  return n;
}

// Constant assignment operation.
template <class nr_type_t>
tvector<nr_type_t> tvector<nr_type_t>::operator = (const nr_type_t val) {
  for (int i = 0; i < size; i++) data[i] = val;
  return *this;
}

// Returns the sum of the vector elements.
template <class nr_type_t>
nr_type_t sum (tvector<nr_type_t> a) {
  nr_type_t res = 0;
  for (int i = 0; i < a.size (); i++) res += a(i);
  return res;
}

// Vector negation.
template <class nr_type_t>
tvector<nr_type_t> operator - (tvector<nr_type_t> a) {
  int n = a.size ();
  tvector<nr_type_t> res (n);
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

// Conjugate vector.
template <class nr_type_t>
tvector<nr_type_t> conj (tvector<nr_type_t> a) {
  int n = a.size ();
  tvector<nr_type_t> res (n);
  for (int i = 0; i < n; i++) res(i) =  conj (a(i));
  return res;
}

// Checks validity of vector.
template <class nr_type_t>
int tvector<nr_type_t>::isFinite (void) {
  for (int i = 0; i < this->size(); i++)
    if (!finite (real (data[i]))) return 0;
  return 1;
}

// The functions reorders the vector according to the given index array.
template <class nr_type_t>
void tvector<nr_type_t>::reorder (int * idx) {
  tvector<nr_type_t> old = *this;
  for (int i = 0; i < this->size(); i++) data[i] = old(idx[i]);
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
