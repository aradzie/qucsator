/*
 * diode.cpp - diode class implementation
 *
 * Copyright (C) 2004 Stefan Jahn <stefan@lkcc.org>
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
 * the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.  
 *
 * $Id: diode.cpp,v 1.6 2004/07/31 16:59:15 ela Exp $
 *
 */

#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "complex.h"
#include "object.h"
#include "node.h"
#include "circuit.h"
#include "net.h"
#include "matrix.h"
#include "analysis.h"
#include "dcsolver.h"
#include "component_id.h"
#include "constants.h"
#include "device.h"
#include "diode.h"

#define NODE_C 1 /* cathode node */
#define NODE_A 2 /* anode node   */

diode::diode () : circuit (2) {
  rs = NULL;
  type = CIR_DIODE;
}

void diode::calcSP (nr_double_t frequency) {
  nr_double_t gd = getOperatingPoint ("gd");
  nr_double_t Cd = getOperatingPoint ("Cd");
  complex y = 2 * z0 * rect (gd, Cd * 2.0 * M_PI * frequency);
  setS (NODE_C, NODE_C, 1.0 / (1.0 + y));
  setS (NODE_A, NODE_A, 1.0 / (1.0 + y));
  setS (NODE_C, NODE_A, y / (1.0 + y));
  setS (NODE_A, NODE_C, y / (1.0 + y));
}

void diode::calcNoise (nr_double_t frequency) {
  nr_double_t Id = getOperatingPoint ("Id");

#if MICHAEL /* shot noise only */
  nr_double_t gd = getOperatingPoint ("gd");
  nr_double_t Cd = getOperatingPoint ("Cd");

  complex y = rect (gd, Cd * 2.0 * M_PI * frequency);
  complex f = 2 * z0 * Id / norm (2 * z0 * y + 1) * QoverkB / T0;
  setN (NODE_C, NODE_C, +f); setN (NODE_A, NODE_A, +f);
  setN (NODE_C, NODE_A, -f); setN (NODE_A, NODE_C, -f);
#else
  nr_double_t Kf  = getPropertyDouble ("Kf");
  nr_double_t Af  = getPropertyDouble ("Af");
  nr_double_t Ffe = getPropertyDouble ("Ffe");

  matrix yc = matrix (2);
  nr_double_t i = 2 * Id * QoverkB / T0 +               // shot noise
    Kf * pow (Id, Af) / pow (frequency, Ffe) / kB / T0; // flicker noise
  yc.set (NODE_C, NODE_C, +i); yc.set (NODE_A, NODE_A, +i);
  yc.set (NODE_A, NODE_C, -i); yc.set (NODE_C, NODE_A, -i);
  setMatrixN (cytocs (yc * z0, getMatrixS ()));
#endif
}

void diode::initDC (dcsolver * solver) {

  // initialize starting values
  setV (NODE_C, 0.0);
  setV (NODE_A, 0.9);
  Uprev = real (getV (NODE_A) - getV (NODE_C));

  // get device temperature
  nr_double_t T = getPropertyDouble ("Temp");

  // possibly insert series resistance
  nr_double_t Rs = getPropertyDouble ("Rs");
  if (Rs != 0) {
    // create additional circuit if necessary and reassign nodes
    rs = splitResistance (this, rs, solver->getNet (), "Rs", "anode", NODE_A);
    rs->setProperty ("Temp", T);
    rs->setProperty ("R", Rs);
  }
  // no series resistance
  else {
    disableResistance (this, rs, solver->getNet (), NODE_A);
  }
}

void diode::calcDC (void) {
  nr_double_t Is = getPropertyDouble ("Is");
  nr_double_t n  = getPropertyDouble ("N");
  nr_double_t T  = getPropertyDouble ("Temp");

  nr_double_t Ud, Ut, Ieq, Ucrit, gtiny;

  T = kelvin (T);
  Ut = T * kBoverQ;
  Ud = real (getV (NODE_A) - getV (NODE_C));

  // critical voltage necessary for bad start values
  Ucrit = pnCriticalVoltage (Is, n * Ut);
  Uprev = Ud = pnVoltage (Ud, Uprev, Ut * n, Ucrit);

  // tiny derivative for little junction voltage
  gtiny = Ud < - 10 * Ut * n ? Is : 0;

  gd = pnConductance (Ud, Is, Ut * n) + gtiny;
  Id = pnCurrent (Ud, Is, Ut * n) + gtiny * Ud;
  Ieq = Id - Ud * gd;

  setI (NODE_C, +Ieq);
  setI (NODE_A, -Ieq);

  setY (NODE_C, NODE_C, +gd); setY (NODE_A, NODE_A, +gd);
  setY (NODE_C, NODE_A, -gd); setY (NODE_A, NODE_C, -gd);
}

void diode::calcOperatingPoints (void) {
  nr_double_t z   = getPropertyDouble ("M");
  nr_double_t cj0 = getPropertyDouble ("Cj0");
  nr_double_t vd  = getPropertyDouble ("Vj");
  nr_double_t Tt  = getPropertyDouble ("Tt");
  
  nr_double_t Ud, Cd;

  Ud = real (getV (NODE_A) - getV (NODE_C));
  Cd = pnCapacitance (Ud, cj0, vd, z) + Tt * gd;

  setOperatingPoint ("gd", gd);
  setOperatingPoint ("Id", Id);
  setOperatingPoint ("Vd", Ud);
  setOperatingPoint ("Cd", Cd);
}
