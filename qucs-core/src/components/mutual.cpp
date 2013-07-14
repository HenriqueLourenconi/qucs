/*
 * mutual.cpp - two mutual inductors class implementation
 *
 * Copyright (C) 2005, 2006, 2008 Stefan Jahn <stefan@lkcc.org>
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
#include "mutual.h"

mutual::mutual () : circuit (4) {
  type = CIR_MUTUAL;
}

void mutual::calcSP (nr_double_t frequency) {
#if 0
  setMatrixS (ytos (calcMatrixY (frequency)));
#else
  nr_double_t l1 = getPropertyDouble ("L1");
  nr_double_t l2 = getPropertyDouble ("L2");
  nr_double_t k = getPropertyDouble ("k");
  nr_double_t o = 2 * M_PI * frequency;
  nr_double_t a = k * k - 1;
  nr_complex_t d = nr_complex_t (o * o * l1 * l2 * a / 2 / z0 + 2 * z0, o * (l1 + l2));
  nr_complex_t r;
  r = nr_complex_t (2 * z0, o * l2) / d;
  setS (NODE_1, NODE_4, r); setS (NODE_4, NODE_1, r);
  r = 1.0 - r;
  setS (NODE_1, NODE_1, r); setS (NODE_4, NODE_4, r);
  r = nr_complex_t (2 * z0, o * l1) / d;
  setS (NODE_2, NODE_3, r); setS (NODE_3, NODE_2, r);
  r = 1.0 - r;
  setS (NODE_2, NODE_2, r); setS (NODE_3, NODE_3, r);
  r = nr_complex_t (0, o * k * sqrt (l1 * l2)) / d;
  setS (NODE_1, NODE_2, r); setS (NODE_2, NODE_1, r);
  setS (NODE_3, NODE_4, r); setS (NODE_4, NODE_3, r);
  r = -r;
  setS (NODE_1, NODE_3, r); setS (NODE_3, NODE_1, r);
  setS (NODE_2, NODE_4, r); setS (NODE_4, NODE_2, r);
#endif
}

matrix mutual::calcMatrixY (nr_double_t frequency) {
  nr_double_t l1 = getPropertyDouble ("L1");
  nr_double_t l2 = getPropertyDouble ("L2");
  nr_double_t k = getPropertyDouble ("k");
  nr_double_t o = 2 * M_PI * frequency;
  nr_double_t a = 1 - k * k;
  nr_complex_t z1 = nr_complex_t (0, o * l1 * a);
  nr_complex_t z2 = nr_complex_t (0, o * l2 * a);
  nr_complex_t y3 = nr_complex_t (0, k / (o * sqrt (l1 * l2) * a));
  
  matrix y = matrix (4);
  y.set (NODE_1, NODE_1, +1.0 / z1); y.set (NODE_4, NODE_4, +1.0 / z1);
  y.set (NODE_1, NODE_4, -1.0 / z1); y.set (NODE_4, NODE_1, -1.0 / z1);
  y.set (NODE_2, NODE_2, +1.0 / z2); y.set (NODE_3, NODE_3, +1.0 / z2);
  y.set (NODE_2, NODE_3, -1.0 / z2); y.set (NODE_3, NODE_2, -1.0 / z2);
  y.set (NODE_1, NODE_3, -y3); y.set (NODE_3, NODE_1, -y3);
  y.set (NODE_2, NODE_4, -y3); y.set (NODE_4, NODE_2, -y3);
  y.set (NODE_1, NODE_2, +y3); y.set (NODE_2, NODE_1, +y3);
  y.set (NODE_3, NODE_4, +y3); y.set (NODE_4, NODE_3, +y3);
  return y;
}

void mutual::initAC (void) {
  setVoltageSources (0);
  allocMatrixMNA ();
}

void mutual::calcAC (nr_double_t frequency) {
  setMatrixY (calcMatrixY (frequency));
}

void mutual::initDC (void) {
  setVoltageSources (2);
  allocMatrixMNA ();
  voltageSource (VSRC_1, NODE_1, NODE_4);
  voltageSource (VSRC_2, NODE_2, NODE_3);
}

void mutual::initTR (void) {
  initDC ();
}

void mutual::calcTR (nr_double_t) {
  nr_double_t k  = getPropertyDouble ("k");
  nr_double_t l1 = getPropertyDouble ("L1");
  nr_double_t l2 = getPropertyDouble ("L2");
  nr_double_t M12 = k * sqrt (l1 * l2);

  setMD (VSRC_1, VSRC_1, -l1); setMD (VSRC_1, VSRC_2, -M12);
  setMD (VSRC_2, VSRC_2, -l2); setMD (VSRC_2, VSRC_1, -M12);
}

// properties
PROP_REQ [] = {
  { "L1", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGEX },
  { "L2", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGEX },
  { "k", PROP_REAL, { 0.9, PROP_NO_STR }, PROP_RNGXX (-1, 1) },
  PROP_NO_PROP };
PROP_OPT [] = {
  PROP_NO_PROP };
struct define_t mutual::cirdef =
  { "MUT", 4, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
