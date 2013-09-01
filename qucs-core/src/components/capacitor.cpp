/*
 * capacitor.cpp - capacitor class implementation
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

/*!\file capacitor.cpp 
   \brief capacitor class implementation

   A capacitor is a passive device that strore electric energy
*/


#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "component.h"
#include "capacitor.h"

/*!\brief Constructor */
capacitor::capacitor () : circuit (2) {
  type = CIR_CAPACITOR;
  setISource (true);
}

/*!\brief Compute S parameters 

  \f$S\f$ parameter are computed from admitance, therefore \f$S\f$
  matrix of a capacitor of capacitance \f$C\f$ is:
  \f[
  S=\begin{pmatrix}
  \frac{1}{1+4j\pi fCZ_0}            & \frac{4j\pi fCZ_0}{1+4j\pi fCZ_0} \\
  \frac{4j\pi fCZ_0}{1+4j\pi fCZ_0}  & \frac{1}{1+4j\pi fCZ_0}
  \end{pmatrix}
  \f]

  \param[in] frequency frequency for S parameters simulation
*/
void capacitor::calcSP (nr_double_t frequency) {
  nr_double_t c = getPropertyDouble ("C") * z0;
  nr_complex_t y = 2.0 * rect (0, 2.0 * M_PI * frequency * c);
  setS (NODE_1, NODE_1, 1.0 / (1.0 + y));
  setS (NODE_2, NODE_2, 1.0 / (1.0 + y));
  setS (NODE_1, NODE_2, y / (1.0 + y));
  setS (NODE_2, NODE_1, y / (1.0 + y));
}

/*\brief Init DC simulation of capacitor */
void capacitor::initDC (void) {
  bool vsource = isPropertyGiven ("V");

  if (vsource)
    setVoltageSources (1);  

  allocMatrixMNA ();

  if (vsource) {
    voltageSource (VSRC_1, NODE_1, NODE_2);
    setE (VSRC_1, getPropertyDouble ("V"));
  }
}

void capacitor::calcDC (void) {
  nr_double_t f = getNet()->getSrcFactor ();
  setE (VSRC_1, f * getPropertyDouble ("V"));
}

/*!\brief AC model 
   
   Capacitor (capacitance \f$C\f$) is modelized by 
   its \f$Y\f$ matrix:
   \f[
   Y=\begin{pmatrix}
     2j\pi f C  & -2j\pi f C \\
     -2j\pi f C & 2j\pi f C
     \end{pmatrix}
   \f]

   \param[in] frequency frequency used for AC simulation
*/
void capacitor::calcAC (nr_double_t frequency) {
  nr_double_t c = getPropertyDouble ("C");
  nr_complex_t y = rect (0, 2.0 * M_PI * frequency * c);
  setY (NODE_1, NODE_1, +y); setY (NODE_2, NODE_2, +y);
  setY (NODE_1, NODE_2, -y); setY (NODE_2, NODE_1, -y);
}

/*!\brief Init AC model of capacitor */
void capacitor::initAC (void) {
  allocMatrixMNA ();
}

void capacitor::initTR (void) {
  setVoltageSources (0);
  allocMatrixMNA ();
}

void capacitor::calcTR (nr_double_t) {

  /* if this is a controlled capacitance then do nothing here */
  if (hasProperty ("Controlled")) return;

  nr_double_t c = getPropertyDouble ("C");

  setMY (NODE_1, NODE_1, +c); setMY (NODE_2, NODE_2, +c);
  setMY (NODE_1, NODE_2, -c); setMY (NODE_2, NODE_1, -c);
}

void capacitor::initHB (void) {
  initAC ();
}

void capacitor::calcHB (nr_double_t frequency) {
  calcAC (frequency);
}

// properties
PROP_REQ [] = {
  { "C", PROP_REAL, { 1e-12, PROP_NO_STR }, PROP_NO_RANGE },
  PROP_NO_PROP };
PROP_OPT [] = {
  { "V", PROP_REAL, { 0, PROP_NO_STR }, PROP_NO_RANGE },
  PROP_NO_PROP };
struct define_t capacitor::cirdef =
  { "C", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
