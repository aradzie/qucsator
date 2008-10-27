/*
 * jfet.cpp - jfet class implementation
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
 * $Id: jfet.cpp,v 1.37 2008-10-27 17:07:17 ela Exp $
 *
 */

#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "component.h"
#include "device.h"
#include "jfet.h"

#define NODE_G 0 /* gate node   */
#define NODE_D 1 /* drain node  */
#define NODE_S 2 /* source node */

using namespace device;

jfet::jfet () : circuit (3) {
  rs = rd = NULL;
  type = CIR_JFET;
}

void jfet::calcSP (nr_double_t frequency) {
  setMatrixS (ytos (calcMatrixY (frequency)));
}

matrix jfet::calcMatrixY (nr_double_t frequency) {

  // fetch computed operating points
  nr_double_t Cgd = getOperatingPoint ("Cgd");
  nr_double_t Cgs = getOperatingPoint ("Cgs");
  nr_double_t ggs = getOperatingPoint ("ggs");
  nr_double_t ggd = getOperatingPoint ("ggd");
  nr_double_t gds = getOperatingPoint ("gds");
  nr_double_t gm  = getOperatingPoint ("gm");

  // compute the models admittances
  nr_complex_t Ygd = rect (ggd, 2.0 * M_PI * frequency * Cgd);
  nr_complex_t Ygs = rect (ggs, 2.0 * M_PI * frequency * Cgs);
  nr_complex_t Yds = gds;

  // build admittance matrix and convert it to S-parameter matrix
  matrix y (3);
  y.set (NODE_G, NODE_G, Ygd + Ygs);
  y.set (NODE_G, NODE_D, -Ygd);
  y.set (NODE_G, NODE_S, -Ygs);
  y.set (NODE_D, NODE_G, gm - Ygd);
  y.set (NODE_D, NODE_D, Ygd + Yds);
  y.set (NODE_D, NODE_S, -Yds - gm);
  y.set (NODE_S, NODE_G, -Ygs - gm);
  y.set (NODE_S, NODE_D, -Yds);
  y.set (NODE_S, NODE_S, Ygs + Yds + gm);
  return y;
}

void jfet::calcNoiseSP (nr_double_t frequency) {
  setMatrixN (cytocs (calcMatrixCy (frequency) * z0, getMatrixS ()));
}

matrix jfet::calcMatrixCy (nr_double_t frequency) {
  /* get operating points and noise properties */
  nr_double_t Kf  = getPropertyDouble ("Kf");
  nr_double_t Af  = getPropertyDouble ("Af");
  nr_double_t Ffe = getPropertyDouble ("Ffe");
  nr_double_t gm  = fabs (getOperatingPoint ("gm"));
  nr_double_t Ids = fabs (getOperatingPoint ("Id"));
  nr_double_t T   = getPropertyDouble ("Temp");

  /* compute channel noise and flicker noise generated by the DC
     transconductance and current flow from drain to source */
  nr_double_t i = 8 * kelvin (T) / T0 * gm / 3 +
    Kf * pow (Ids, Af) / pow (frequency, Ffe) / kB / T0;

  /* build noise current correlation matrix and convert it to
     noise-wave correlation matrix */
  matrix cy = matrix (3);
  cy.set (NODE_D, NODE_D, +i);
  cy.set (NODE_S, NODE_S, +i);
  cy.set (NODE_D, NODE_S, -i);
  cy.set (NODE_S, NODE_D, -i);
  return cy;
}

void jfet::initModel (void) {
  // fetch necessary device properties
  nr_double_t T  = getPropertyDouble ("Temp");
  nr_double_t Tn = getPropertyDouble ("Tnom");
  nr_double_t A  = getPropertyDouble ("Area");

  // compute Is temperature and area dependency
  nr_double_t Is  = getPropertyDouble ("Is");
  nr_double_t N   = getPropertyDouble ("N");
  nr_double_t Xti = getPropertyDouble ("Xti");
  nr_double_t T1, T2, Eg;
  T2 = kelvin (T);
  T1 = kelvin (Tn);
  Eg = Egap (300);
  Is = pnCurrent_T (T1, T2, Is, Eg, N, Xti);
  setScaledProperty ("Is", Is * A);

  // compute Isr temperature and area dependency
  nr_double_t Isr = getPropertyDouble ("Isr");
  nr_double_t Nr  = getPropertyDouble ("Nr");
  Isr = pnCurrent_T (T1, T2, Isr, Eg, Nr, Xti);
  setScaledProperty ("Isr", Isr * A);

  // compute Pb temperature dependency
  nr_double_t Pb = getPropertyDouble ("Pb");
  nr_double_t PbT;
  PbT = pnPotential_T (T1,T2, Pb);
  setScaledProperty ("Pb", PbT);

  // compute Cgs and Cgd temperature and area dependency
  nr_double_t Cgs = getPropertyDouble ("Cgs");
  nr_double_t Cgd = getPropertyDouble ("Cgd");
  nr_double_t M   = getPropertyDouble ("M");
  nr_double_t F;
  F = A * pnCapacitance_F (T1, T2, M, PbT / Pb);
  setScaledProperty ("Cgs", Cgs * F);
  setScaledProperty ("Cgd", Cgd * F);

  // compute Vth temperature dependency
  nr_double_t Vt0   = getPropertyDouble ("Vt0");
  nr_double_t Vt0tc = getPropertyDouble ("Vt0tc");
  nr_double_t DT    = T2 - T1;
  Vt0 = Vt0 + Vt0tc * DT;
  setScaledProperty ("Vt0", Vt0);

  // compute Beta temperature and area dependency
  nr_double_t Beta    = getPropertyDouble ("Beta");
  nr_double_t Betatce = getPropertyDouble ("Betatce");
  Beta = Beta * exp (Betatce * DT * log (1.01));
  setScaledProperty ("Beta", Beta * A);

  // compute Rs and Rd area dependency
  nr_double_t Rs = getPropertyDouble ("Rs");
  nr_double_t Rd = getPropertyDouble ("Rd");
  setScaledProperty ("Rs", Rs / A);
  setScaledProperty ("Rd", Rd / A);
}

void jfet::restartDC (void) {
  // apply starting values to previous iteration values
  UgdPrev = real (getV (NODE_G) - getV (NODE_D));
  UgsPrev = real (getV (NODE_G) - getV (NODE_S));
}

void jfet::initDC (void) {

  // allocate MNA matrices
  allocMatrixMNA ();

  // initialize scalability
  initModel ();

  // initialize starting values
  restartDC ();

  // apply polarity of JFET
  char * type = getPropertyString ("Type");
  pol = !strcmp (type, "pfet") ? -1 : 1;

  // get device temperature
  nr_double_t T = getPropertyDouble ("Temp");

  // possibly insert series resistance at source
  nr_double_t Rs = getScaledProperty ("Rs");
  if (Rs != 0.0) {
    // create additional circuit if necessary and reassign nodes
    rs = splitResistor (this, rs, "Rs", "source", NODE_S);
    rs->setProperty ("Temp", T);
    rs->setProperty ("R", Rs);
    rs->setProperty ("Controlled", getName ());
    rs->initDC ();
  }
  // no series resistance at source
  else {
    disableResistor (this, rs, NODE_S);
  }

  // possibly insert series resistance at drain
  nr_double_t Rd = getScaledProperty ("Rd");
  if (Rd != 0.0) {
    // create additional circuit if necessary and reassign nodes
    rd = splitResistor (this, rd, "Rd", "drain", NODE_D);
    rd->setProperty ("Temp", T);
    rd->setProperty ("R", Rd);
    rd->setProperty ("Controlled", getName ());
    rd->initDC ();
  }
  // no series resistance at drain
  else {
    disableResistor (this, rd, NODE_D);
  }
}

void jfet::calcDC (void) {

  // fetch device model parameters
  nr_double_t Is   = getScaledProperty ("Is");
  nr_double_t n    = getPropertyDouble ("N");
  nr_double_t Isr  = getScaledProperty ("Isr");
  nr_double_t nr   = getPropertyDouble ("Nr");
  nr_double_t Vt0  = getScaledProperty ("Vt0");
  nr_double_t l    = getPropertyDouble ("Lambda");
  nr_double_t beta = getScaledProperty ("Beta");
  nr_double_t T    = getPropertyDouble ("Temp");

  nr_double_t Ut, IeqG, IeqD, IeqS, UgsCrit, UgdCrit;
  nr_double_t Igs, Igd, gtiny;

  T = kelvin (T);
  Ut = T * kBoverQ;
  Ugd = real (getV (NODE_G) - getV (NODE_D)) * pol;
  Ugs = real (getV (NODE_G) - getV (NODE_S)) * pol;

  // critical voltage necessary for bad start values
  UgsCrit = pnCriticalVoltage (Is, Ut * n);
  UgdCrit = pnCriticalVoltage (Is, Ut * n);
  UgsPrev = Ugs = pnVoltage (Ugs, UgsPrev, Ut * n, UgsCrit);
  UgdPrev = Ugd = pnVoltage (Ugd, UgdPrev, Ut * n, UgdCrit);

  Uds = Ugs - Ugd;

  // gate-source diode
  gtiny = Ugs < - 10 * Ut * n ? (Is + Isr) : 0;
  ggs = pnConductance (Ugs, Is, Ut * n) +
    pnConductance (Ugs, Isr, Ut * nr) + gtiny;
  Igs = pnCurrent (Ugs, Is, Ut * n) +
    pnCurrent (Ugs, Isr, Ut * nr) + gtiny * Ugs;

  // gate-drain diode
  gtiny = Ugd < - 10 * Ut * n ? (Is + Isr) : 0;
  ggd = pnConductance (Ugd, Is, Ut * n) +
    pnConductance (Ugd, Isr, Ut * nr) + gtiny;
  Igd = pnCurrent (Ugd, Is, Ut * n) +
    pnCurrent (Ugd, Isr, Ut * nr) + gtiny * Ugd;

  // normal (forward) mode of operation
  if (Uds >= 0) {
    nr_double_t Ugst = Ugs - Vt0;
    // normal mode, cutoff region
    if (Ugst <= 0) {
      Ids = 0;
      gm  = 0;
      gds = 0;
    }
    else {
      nr_double_t b = beta * (1 + l * Uds);
      // normal mode, saturation region
      if (Ugst <= Uds) {
	Ids = b * Ugst * Ugst;
	gm  = b * 2 * Ugst;
	gds = l * beta * Ugst * Ugst;
      }
      // normal mode, linear region
      else {
	Ids = b * Uds * (2 * Ugst - Uds);
	gm  = b * 2 * Uds;
	gds = b * 2 * (Ugst - Uds) + l * beta * Uds * (2 * Ugst - Uds);
      }
    }
  }
  // inverse (backward) mode of operation
  else {
    nr_double_t Ugdt = Ugd - Vt0;
    // inverse mode, cutoff region
    if (Ugdt <= 0) {
      Ids = 0;
      gm  = 0;
      gds = 0;
    }
    else {
      nr_double_t b = beta * (1 - l * Uds);
      // inverse mode, saturation region
      if (Ugdt <= -Uds) {
	Ids = - b * Ugdt * Ugdt;
	gm  = - b * 2 * Ugdt;
	gds = beta * l * Ugdt * Ugdt + b * 2 * Ugdt;
      }
      // inverse mode, linear region
      else {
	Ids = b * Uds * (2 * Ugdt + Uds);
	gm  = b * 2 * Uds;
	gds = 2 * b * Ugdt - beta * l * Uds * (2 * Ugdt + Uds);
      }
    }
  }

  // compute autonomic current sources
  IeqG = Igs - ggs * Ugs;
  IeqD = Igd - ggd * Ugd;
  IeqS = Ids - gm * Ugs - gds * Uds;
  setI (NODE_G, (-IeqG - IeqD) * pol);
  setI (NODE_D, (+IeqD - IeqS) * pol);
  setI (NODE_S, (+IeqG + IeqS) * pol);

  // apply admittance matrix elements
  setY (NODE_G, NODE_G, ggs + ggd);
  setY (NODE_G, NODE_D, -ggd);
  setY (NODE_G, NODE_S, -ggs);
  setY (NODE_D, NODE_G, -ggd + gm);
  setY (NODE_D, NODE_D, gds + ggd);
  setY (NODE_D, NODE_S, -gm - gds);
  setY (NODE_S, NODE_G, -ggs - gm);
  setY (NODE_S, NODE_D, -gds);
  setY (NODE_S, NODE_S, ggs + gds + gm);
}

void jfet::loadOperatingPoints (void) {
  Ugs = getOperatingPoint ("Vgs");
  Ugd = getOperatingPoint ("Vgd");
  Uds = getOperatingPoint ("Vds");
}

void jfet::saveOperatingPoints (void) {
  nr_double_t Vgs, Vgd;
  Vgd = real (getV (NODE_G) - getV (NODE_D)) * pol;
  Vgs = real (getV (NODE_G) - getV (NODE_S)) * pol;
  setOperatingPoint ("Vgs", Vgs);
  setOperatingPoint ("Vgd", Vgd);
  setOperatingPoint ("Vds", Vgs - Vgd);
}

void jfet::calcOperatingPoints (void) {

  // fetch device model parameters
  nr_double_t z    = getPropertyDouble ("M");
  nr_double_t Cgd0 = getScaledProperty ("Cgd");
  nr_double_t Cgs0 = getScaledProperty ("Cgs");
  nr_double_t Pb   = getScaledProperty ("Pb");
  nr_double_t Fc   = getPropertyDouble ("Fc");
  
  nr_double_t Cgs, Cgd;

  // capacitance of gate-drain diode
  Cgd = pnCapacitance (Ugd, Cgd0, Pb, z, Fc);
  Qgd = pnCharge (Ugd, Cgd0, Pb, z, Fc);

  // capacitance of gate-source diode
  Cgs = pnCapacitance (Ugs, Cgs0, Pb, z, Fc);
  Qgs = pnCharge (Ugs, Cgs0, Pb, z, Fc);

  // save operating points
  setOperatingPoint ("ggs", ggs);
  setOperatingPoint ("ggd", ggd);
  setOperatingPoint ("gds", gds);
  setOperatingPoint ("gm", gm);
  setOperatingPoint ("Id", Ids);
  setOperatingPoint ("Cgd", Cgd);
  setOperatingPoint ("Cgs", Cgs);
}

void jfet::initAC (void) {
  allocMatrixMNA ();
  clearI ();
}

void jfet::calcAC (nr_double_t frequency) {
  setMatrixY (calcMatrixY (frequency));
}

void jfet::calcNoiseAC (nr_double_t frequency) {
  setMatrixN (calcMatrixCy (frequency));
}

#define qgdState 0 // gate-drain charge state
#define cgdState 1 // gate-drain current state
#define qgsState 2 // gate-source charge state
#define cgsState 3 // gate-source current state

void jfet::initTR (void) {
  setStates (4);
  initDC ();
}

void jfet::calcTR (nr_double_t) {
  calcDC ();
  saveOperatingPoints ();
  loadOperatingPoints ();
  calcOperatingPoints ();

  nr_double_t Cgs = getOperatingPoint ("Cgs");
  nr_double_t Cgd = getOperatingPoint ("Cgd");

  transientCapacitance (qgsState, NODE_G, NODE_S, Cgs, Ugs, Qgs);
  transientCapacitance (qgdState, NODE_G, NODE_D, Cgd, Ugd, Qgd);
}

// properties
PROP_REQ [] = {
  { "Is", PROP_REAL, { 1e-14, PROP_NO_STR }, PROP_POS_RANGE },
  { "N", PROP_REAL, { 1, PROP_NO_STR }, PROP_RNGII (1, 100) },
  { "Vt0", PROP_REAL, { -2, PROP_NO_STR }, PROP_NEG_RANGE },
  { "Lambda", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Beta", PROP_REAL, { 1e-4, PROP_NO_STR }, PROP_POS_RANGE },
  { "M", PROP_REAL, { 0.5, PROP_NO_STR }, PROP_RNGII (0, 1) },
  { "Pb", PROP_REAL, { 1.0, PROP_NO_STR }, PROP_RNGXI (0, 10) },
  { "Fc", PROP_REAL, { 0.5, PROP_NO_STR }, PROP_RNGIX (0, 1) },
  { "Cgs", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Cgd", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  PROP_NO_PROP };
PROP_OPT [] = {
  { "Rd", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Rs", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Isr", PROP_REAL, { 1e-14, PROP_NO_STR }, PROP_POS_RANGE },
  { "Nr", PROP_REAL, { 2, PROP_NO_STR }, PROP_RNGII (1, 100) },
  { "Kf", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Af", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGE },
  { "Ffe", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGE },
  { "Temp", PROP_REAL, { 26.85, PROP_NO_STR }, PROP_MIN_VAL (K) },
  { "Type", PROP_STR, { PROP_NO_VAL, "nfet" }, PROP_RNG_FET },
  { "Xti", PROP_REAL, { 3, PROP_NO_STR }, PROP_POS_RANGE },
  { "Vt0tc", PROP_REAL, { 0, PROP_NO_STR }, PROP_NO_RANGE },
  { "Betatce", PROP_REAL, { 0, PROP_NO_STR }, PROP_NO_RANGE },
  { "Tnom", PROP_REAL, { 26.85, PROP_NO_STR }, PROP_MIN_VAL (K) },
  { "Area", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGEX },
  PROP_NO_PROP };
struct define_t jfet::cirdef =
  { "JFET", 3, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_NONLINEAR, PROP_DEF };
