/*
 * mutual2.cpp - three mutual inductors class implementation
 *
 * Copyright (C) 2005, 2008 Stefan Jahn <stefan@lkcc.org>
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

#include "component.h"
#include "mutual2.h"

mutual2::mutual2 () : circuit (6) {
  type = CIR_MUTUAL2;
}

void mutual2::calcSP (nr_double_t frequency) {
  setMatrixS (ytos (calcMatrixY (frequency)));
}

matrix mutual2::calcMatrixY (nr_double_t frequency) {
  nr_double_t k12 = getPropertyDouble ("k12");
  nr_double_t k13 = getPropertyDouble ("k13");
  nr_double_t k23 = getPropertyDouble ("k23");
  nr_double_t l1 = getPropertyDouble ("L1");
  nr_double_t l2 = getPropertyDouble ("L2");
  nr_double_t l3 = getPropertyDouble ("L3");
  nr_double_t o = 2 * M_PI * frequency;
  nr_double_t a = 1 - k12 * k12 - k13 * k13 - k23 * k23 + 2 * k12 * k13 * k23;
  nr_complex_t y11 = nr_complex_t (0, (k23 * k23 - 1) / l1 / a / o);
  nr_complex_t y22 = nr_complex_t (0, (k12 * k12 - 1) / l3 / a / o);
  nr_complex_t y44 = nr_complex_t (0, (k13 * k13 - 1) / l2 / a / o);
  nr_complex_t y12 = nr_complex_t (0, (k13 - k12 * k23) / sqrt (l1 * l3) / a / o);
  nr_complex_t y15 = nr_complex_t (0, (k12 - k13 * k23) / sqrt (l1 * l2) / a / o);
  nr_complex_t y25 = nr_complex_t (0, (k23 - k12 * k13) / sqrt (l2 * l3) / a / o);

  matrix y = matrix (6);
  y.set (NODE_1, NODE_1, +y11); y.set (NODE_6, NODE_6, +y11);
  y.set (NODE_1, NODE_6, -y11); y.set (NODE_6, NODE_1, -y11);
  y.set (NODE_2, NODE_2, +y22); y.set (NODE_3, NODE_3, +y22);
  y.set (NODE_2, NODE_3, -y22); y.set (NODE_3, NODE_2, -y22);
  y.set (NODE_4, NODE_4, +y44); y.set (NODE_5, NODE_5, +y44);
  y.set (NODE_4, NODE_5, -y44); y.set (NODE_5, NODE_4, -y44);
  y.set (NODE_1, NODE_2, +y12); y.set (NODE_2, NODE_1, +y12);
  y.set (NODE_3, NODE_6, +y12); y.set (NODE_6, NODE_3, +y12);
  y.set (NODE_1, NODE_3, -y12); y.set (NODE_3, NODE_1, -y12);
  y.set (NODE_2, NODE_6, -y12); y.set (NODE_6, NODE_2, -y12);
  y.set (NODE_1, NODE_5, +y15); y.set (NODE_5, NODE_1, +y15);
  y.set (NODE_4, NODE_6, +y15); y.set (NODE_6, NODE_4, +y15);
  y.set (NODE_1, NODE_4, -y15); y.set (NODE_4, NODE_1, -y15);
  y.set (NODE_5, NODE_6, -y15); y.set (NODE_6, NODE_5, -y15);
  y.set (NODE_2, NODE_5, +y25); y.set (NODE_5, NODE_2, +y25);
  y.set (NODE_4, NODE_3, +y25); y.set (NODE_3, NODE_4, +y25);
  y.set (NODE_2, NODE_4, -y25); y.set (NODE_4, NODE_2, -y25);
  y.set (NODE_5, NODE_3, -y25); y.set (NODE_3, NODE_5, -y25);
  return y;
}

void mutual2::initAC (void) {
  setVoltageSources (0);
  allocMatrixMNA ();
}

void mutual2::calcAC (nr_double_t frequency) {
  setMatrixY (calcMatrixY (frequency));
}

void mutual2::initDC (void) {
  setVoltageSources (3);
  allocMatrixMNA ();
  voltageSource (VSRC_1, NODE_1, NODE_6);
  voltageSource (VSRC_2, NODE_5, NODE_4);
  voltageSource (VSRC_3, NODE_2, NODE_3);
}

void mutual2::initTR (void) {
  initDC ();
  setStates (18);
}

void mutual2::calcTR (nr_double_t) {
  nr_double_t k12 = getPropertyDouble ("k12");
  nr_double_t k13 = getPropertyDouble ("k13");
  nr_double_t k23 = getPropertyDouble ("k23");
  nr_double_t l1  = getPropertyDouble ("L1");
  nr_double_t l2  = getPropertyDouble ("L2");
  nr_double_t l3  = getPropertyDouble ("L3");
  nr_double_t M12 = k12 * sqrt (l1 * l2);
  nr_double_t M13 = k13 * sqrt (l1 * l3);
  nr_double_t M23 = k23 * sqrt (l2 * l3);

  setMD (VSRC_1, VSRC_1, -l1);
  setMD (VSRC_1, VSRC_2, -M12);
  setMD (VSRC_1, VSRC_3, -M13);
  setMD (VSRC_2, VSRC_1, -M12);
  setMD (VSRC_2, VSRC_2, -l2);
  setMD (VSRC_2, VSRC_3, -M23);
  setMD (VSRC_3, VSRC_1, -M13);
  setMD (VSRC_3, VSRC_2, -M23);
  setMD (VSRC_3, VSRC_3, -l3);
}

// properties
PROP_REQ [] = {
  { "L1", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGEX },
  { "L2", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGEX },
  { "L3", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGEX },
  { "k12", PROP_REAL, { 0.9, PROP_NO_STR }, PROP_RNGXX (-1, 1) },
  { "k13", PROP_REAL, { 0.9, PROP_NO_STR }, PROP_RNGXX (-1, 1) },
  { "k23", PROP_REAL, { 0.9, PROP_NO_STR }, PROP_RNGXX (-1, 1) },
  PROP_NO_PROP };
PROP_OPT [] = {
  PROP_NO_PROP };
struct define_t mutual2::cirdef =
  { "MUT2", 6, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
