/*
 * inductor.cpp - inductor class implementation
 *
 * Copyright (C) 2003, 2004, 2005, 2006, 2008 Stefan Jahn <stefan@lkcc.org>
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
#include "inductor.h"

inductor::inductor () : circuit (2) {
  type = CIR_INDUCTOR;
  setISource (true);
}

void inductor::calcSP (nr_double_t frequency) {
  nr_double_t l = getPropertyDouble ("L") / z0;
  nr_complex_t z = nr_complex_t (0, 2.0 * M_PI * frequency * l);
  setS (NODE_1, NODE_1, z / (z + 2.0));
  setS (NODE_2, NODE_2, z / (z + 2.0));
  setS (NODE_1, NODE_2, 2.0 / (z + 2.0));
  setS (NODE_2, NODE_1, 2.0 / (z + 2.0));
}

void inductor::initDC (void) {
  if (isPropertyGiven ("I"))
    initDC_isrc ();
  else
    initDC_short ();
}

void inductor::initDC_short (void) {
  setVoltageSources (1);  
  allocMatrixMNA ();
  voltageSource (VSRC_1, NODE_1, NODE_2);
}

void inductor::initDC_isrc (void) {
  nr_double_t i = getPropertyDouble ("I");
  setVoltageSources (1);
  allocMatrixMNA ();
  setI (NODE_1, -i); setI (NODE_2, +i);
  setD (VSRC_1, VSRC_1, 1.0);
  setE (VSRC_1, i);
}

void inductor::calcDC (void) {
  clearY ();
  if (isPropertyGiven ("I")) {
    nr_double_t i = getPropertyDouble ("I");
    nr_double_t f = getNet()->getSrcFactor ();
    setI (NODE_1, -i * f); setI (NODE_2, +i * f);
  }
}

void inductor::initAC (void) {
  nr_double_t l = getPropertyDouble ("L");

  // for non-zero inductance usual MNA entries
  if (l != 0.0) {
    setVoltageSources (0);
    allocMatrixMNA ();
  }
  // for zero inductance create a zero voltage source
  else {
    initDC_short ();
    clearY ();
  }
}

void inductor::calcAC (nr_double_t frequency) {
  nr_double_t l = getPropertyDouble ("L");

  // for non-zero inductance usual MNA entries
  if (l != 0.0) {
    nr_complex_t y = nr_complex_t (0, -1 / (2.0 * M_PI * frequency * l));
    setY (NODE_1, NODE_1, +y); setY (NODE_2, NODE_2, +y);
    setY (NODE_1, NODE_2, -y); setY (NODE_2, NODE_1, -y);
  }
}

void inductor::initTR (void) {
  initDC_short ();
  clearY ();
}

void inductor::calcTR (nr_double_t) {
  nr_double_t l = getPropertyDouble ("L");

  setMD (VSRC_1, VSRC_1, -l);
}

void inductor::initHB (void) {
  setVoltageSources (1);
  setInternalVoltageSource (1);
  allocMatrixMNA ();
  voltageSource (VSRC_1, NODE_1, NODE_2);
}

void inductor::calcHB (nr_double_t frequency) {
  nr_double_t l = getPropertyDouble ("L");
  setD (VSRC_1, VSRC_1, -l * 2 * M_PI * frequency);
}

// properties
PROP_REQ [] = {
  { "L", PROP_REAL, { 1e-9, PROP_NO_STR }, PROP_NO_RANGE }, PROP_NO_PROP };
PROP_OPT [] = {
  { "I", PROP_REAL, { 0, PROP_NO_STR }, PROP_NO_RANGE }, PROP_NO_PROP };
struct define_t inductor::cirdef =
  { "L", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
