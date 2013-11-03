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

template <typename nr_type_t>
using tmatrix = Eigen::Matrix<nr_type_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;


#endif /* __TMATRIX_H__ */
