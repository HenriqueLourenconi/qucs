/*
 * tmatrix.cpp - simple matrix template class implementation
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
#include "logging.h"
#include "complex.h"
#include "tmatrix.h"

#ifdef DEBUG
// Debug function: Prints the matrix object.
template <class nr_type_t>
void tmatrix<nr_type_t>::print (bool realonly) {
  for (int r = 0; r < this->m.rows; r++) {
    for (int c = 0; c < this->m.cols(); c++) {
      if (realonly) {
	fprintf (stderr, "%+.2e%s", (double) real (this->m(r, c)),
		 c != this->m.cols() - 1 ? " " : "");
      } else {
	fprintf (stderr, "%+.2e%+.2ei%s", (double) real (this->m(r, c)),
		 (double) imag (this->m(r, c)), c !=  this->m.cols() - 1 ? " " : "");
      }
    }
    fprintf (stderr, ";\n");
  }
  fprintf (stderr, ";;;\n");
  fflush (stderr);
}
#endif /* DEBUG */

