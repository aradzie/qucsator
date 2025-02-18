/*
 * bjt.cpp - bipolar junction transistor class implementation
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
 */

#include "component.h"
#include "device.h"
#include "bjt.h"

#define NEWSGP 0

#define NODE_B 0 /* base node       */
#define NODE_C 1 /* collector node  */
#define NODE_E 2 /* emitter node    */
#define NODE_S 3 /* substrate node  */

using namespace qucs;
using namespace qucs::device;

bjt::bjt () : circuit (4) {
  cbcx = rb = re = rc = NULL;
  type = CIR_BJT;
}

void bjt::calcSP (double frequency) {
  // build admittance matrix and convert it to S-parameter matrix
  setMatrixS (ytos (calcMatrixY (frequency)));
}

matrix bjt::calcMatrixY (double frequency) {

  // fetch computed operating points
  double Cbe  = getOperatingPoint ("Cbe");
  double gbe  = getOperatingPoint ("gpi");
  double Cbci = getOperatingPoint ("Cbci");
  double gbc  = getOperatingPoint ("gmu");
  double Ccs  = getOperatingPoint ("Ccs");
#if NEWSGP
  double gm   = getOperatingPoint ("gmf");
  double gmr  = getOperatingPoint ("gmr");
#else
  double gm   = getOperatingPoint ("gm");
  double go   = getOperatingPoint ("go");
#endif
  double Ptf  = getPropertyDouble ("Ptf");
  double Tf   = getPropertyDouble ("Tf");

  // compute admittance matrix entries
  nr_complex_t Ybe = nr_complex_t (gbe, 2.0 * pi * frequency * Cbe);
  nr_complex_t Ybc = nr_complex_t (gbc, 2.0 * pi * frequency * Cbci);
  nr_complex_t Ycs = nr_complex_t (0.0, 2.0 * pi * frequency * Ccs);

  // admittance matrix entries for "transcapacitance"
  nr_complex_t Ybebc = nr_complex_t (0.0, 2.0 * pi * frequency * dQbedUbc);

  // compute influence of excess phase
  double phase = deg2rad (Ptf) * Tf * 2 * pi * frequency;
#if NEWSGP
  nr_complex_t gmf = qucs::polar (gm, -phase);
#else
  nr_complex_t gmf = qucs::polar (gm + go, -phase) - go;
#endif

  // build admittance matrix
  matrix y (4);
#if NEWSGP
  // for some reason this small signal equivalent can't be used
  y.set (NODE_B, NODE_B, Ybc + Ybe + Ybebc);
  y.set (NODE_B, NODE_C, -Ybc - Ybebc);
  y.set (NODE_B, NODE_E, -Ybe);
  y.set (NODE_B, NODE_S, 0);
  y.set (NODE_C, NODE_B, -Ybc + gmf + gmr);
  y.set (NODE_C, NODE_C, Ybc - gmr + Ycs);
  y.set (NODE_C, NODE_E, -gmf);
  y.set (NODE_C, NODE_S, -Ycs);
  y.set (NODE_E, NODE_B, -Ybe - gmf - gmr - Ybebc);
  y.set (NODE_E, NODE_C, gmr + Ybebc);
  y.set (NODE_E, NODE_E, Ybe + gmf);
  y.set (NODE_E, NODE_S, 0);
  y.set (NODE_S, NODE_B, 0);
  y.set (NODE_S, NODE_C, -Ycs);
  y.set (NODE_S, NODE_E, 0);
  y.set (NODE_S, NODE_S, Ycs);
#else /* !NEWSGP */
  y.set (NODE_B, NODE_B, Ybc + Ybe + Ybebc);
  y.set (NODE_B, NODE_C, -Ybc - Ybebc);
  y.set (NODE_B, NODE_E, -Ybe);
  y.set (NODE_B, NODE_S, 0);
  y.set (NODE_C, NODE_B, -Ybc + gmf);
  y.set (NODE_C, NODE_C, Ybc + Ycs + go);
  y.set (NODE_C, NODE_E, -gmf - go);
  y.set (NODE_C, NODE_S, -Ycs);
  y.set (NODE_E, NODE_B, -Ybe - gmf - Ybebc);
  y.set (NODE_E, NODE_C, -go + Ybebc);
  y.set (NODE_E, NODE_E, Ybe + gmf + go);
  y.set (NODE_E, NODE_S, 0);
  y.set (NODE_S, NODE_B, 0);
  y.set (NODE_S, NODE_C, -Ycs);
  y.set (NODE_S, NODE_E, 0);
  y.set (NODE_S, NODE_S, Ycs);
#endif /* !NEWSGP */
  return y;
}

void bjt::calcNoiseSP (double frequency) {
  setMatrixN (cytocs (calcMatrixCy (frequency) * z0, getMatrixS ()));
}

matrix bjt::calcMatrixCy (double frequency) {

  // fetch computed operating points
  double Ibe = fabs (getOperatingPoint ("Ibe"));
  double Ice = fabs (getOperatingPoint ("Ice"));

  // get model properties
  double Kf  = getPropertyDouble ("Kf");
  double Af  = getPropertyDouble ("Af");
  double Ffe = getPropertyDouble ("Ffe");
  double Kb  = getPropertyDouble ("Kb");
  double Ab  = getPropertyDouble ("Ab");
  double Fb  = getPropertyDouble ("Fb");

  double ib = 2 * Ibe * QoverkB / T0 +            // shot noise
    (Kf * qucs::pow (Ibe, Af) / qucs::pow (frequency, Ffe) +       // flicker noise
     Kb * qucs::pow (Ibe, Ab) / (1 + sqr (frequency / Fb)))  // burst noise
    / kB / T0;
  double ic = 2 * Ice * QoverkB / T0;             // shot noise

  /* build noise current correlation matrix and convert it to
     noise-wave correlation matrix */
  matrix cy = matrix (4);
  cy.set (NODE_B, NODE_B, ib);
  cy.set (NODE_B, NODE_E, -ib);
  cy.set (NODE_C, NODE_C, ic);
  cy.set (NODE_C, NODE_E, -ic);
  cy.set (NODE_E, NODE_B, -ib);
  cy.set (NODE_E, NODE_C, -ic);
  cy.set (NODE_E, NODE_E, ic + ib);
  return cy;
}

void bjt::initModel (void) {
  // fetch necessary device properties
  double T  = getPropertyDouble ("Temp");
  double Tn = getPropertyDouble ("Tnom");
  double A  = getPropertyDouble ("Area");

  // compute Is temperature and area dependency
  double Is  = getPropertyDouble ("Is");
  double Xti = getPropertyDouble ("Xti");
  double Eg  = getPropertyDouble ("Eg");
  double T1, T2, IsT;
  T2 = celsius2kelvin (T);
  T1 = celsius2kelvin (Tn);
  IsT = pnCurrent_T (T1, T2, Is, Eg, 1, Xti);
  setScaledProperty ("Is", IsT * A);

  // compute Vje, Vjc and Vjs temperature dependencies
  double Vje = getPropertyDouble ("Vje");
  double Vjc = getPropertyDouble ("Vjc");
  double Vjs = getPropertyDouble ("Vjs");
  double VjeT, VjcT, VjsT;
  VjeT = pnPotential_T (T1,T2, Vje);
  VjcT = pnPotential_T (T1,T2, Vjc);
  VjsT = pnPotential_T (T1,T2, Vjs);
  setScaledProperty ("Vje", VjeT);
  setScaledProperty ("Vjc", VjcT);
  setScaledProperty ("Vjs", VjsT);

  // compute Bf and Br temperature dependencies
  double Bf  = getPropertyDouble ("Bf");
  double Br  = getPropertyDouble ("Br");
  double Xtb = getPropertyDouble ("Xtb");
  double F = qucs::exp (Xtb * qucs::log (T2 / T1));
  setScaledProperty ("Bf", Bf * F);
  setScaledProperty ("Br", Br * F);

  // compute Ise and Isc temperature and area dependencies
  double Ise = getPropertyDouble ("Ise");
  double Isc = getPropertyDouble ("Isc");
  double Ne  = getPropertyDouble ("Ne");
  double Nc  = getPropertyDouble ("Nc");
  double G = qucs::log (IsT / Is);
  double F1 = qucs::exp (G / Ne);
  double F2 = qucs::exp (G / Nc);
  Ise = Ise / F * F1;
  Isc = Isc / F * F2;
  setScaledProperty ("Ise", Ise * A);
  setScaledProperty ("Isc", Isc * A);

  // check unphysical parameters
  double Nf = getPropertyDouble ("Nf");
  double Nr = getPropertyDouble ("Nr");
  if (Nf < 1.0) {
    logprint (LOG_ERROR, "WARNING: Unphysical model parameter Nf = %g in "
	      "BJT `%s'\n", Nf, getName ());
  }
  if (Nr < 1.0) {
    logprint (LOG_ERROR, "WARNING: Unphysical model parameter Nr = %g in "
	      "BJT `%s'\n", Nr, getName ());
  }
  if (Ne < 1.0) {
    logprint (LOG_ERROR, "WARNING: Unphysical model parameter Ne = %g in "
	      "BJT `%s'\n", Ne, getName ());
  }
  if (Nc < 1.0) {
    logprint (LOG_ERROR, "WARNING: Unphysical model parameter Nc = %g in "
	      "BJT `%s'\n", Nc, getName ());
  }

  /* Originally Vtf was expected to be PROP_POS_RANGE, but there are models
   * which use negative values. Instead of failing, warn the user.
   * \todo Provide a way to silece such warnings
   */
  double Vtf = getPropertyDouble ("Vtf");
  if (Vtf < 0.0) {
    logprint (LOG_ERROR, "WARNING: Unphysical model parameter Vtf = %g in "
	      "BJT `%s'\n", Vtf, getName ());
  }

  // compute Cje, Cjc and Cjs temperature and area dependencies
  double Cje = getPropertyDouble ("Cje");
  double Cjc = getPropertyDouble ("Cjc");
  double Cjs = getPropertyDouble ("Cjs");
  double Mje = getPropertyDouble ("Mje");
  double Mjc = getPropertyDouble ("Mjc");
  double Mjs = getPropertyDouble ("Mjs");
  Cje = pnCapacitance_T (T1, T2, Mje, VjeT / Vje, Cje);
  Cjc = pnCapacitance_T (T1, T2, Mjc, VjcT / Vjc, Cjc);
  Cjs = pnCapacitance_T (T1, T2, Mjs, VjsT / Vjs, Cjs);
  setScaledProperty ("Cje", Cje * A);
  setScaledProperty ("Cjc", Cjc * A);
  setScaledProperty ("Cjs", Cjs * A);

  // compute Rb, Rc, Re and Rbm area dependencies
  double Rb  = getPropertyDouble ("Rb");
  double Re  = getPropertyDouble ("Re");
  double Rc  = getPropertyDouble ("Rc");
  double Rbm = getPropertyDouble ("Rbm");
  setScaledProperty ("Rb", Rb / A);
  setScaledProperty ("Re", Re / A);
  setScaledProperty ("Rc", Rc / A);
  setScaledProperty ("Rbm", Rbm / A);

  // compute Ikf, Ikr, Irb and Itf area dependencies
  double Ikf = getPropertyDouble ("Ikf");
  double Ikr = getPropertyDouble ("Ikr");
  double Irb = getPropertyDouble ("Irb");
  double Itf = getPropertyDouble ("Itf");
  setScaledProperty ("Ikf", Ikf * A);
  setScaledProperty ("Ikr", Ikr * A);
  setScaledProperty ("Irb", Irb * A);
  setScaledProperty ("Itf", Itf * A);
}

void bjt::initDC (void) {

  // no transient analysis
  doTR = false;

  // allocate MNA matrices
  allocMatrixMNA ();

  // initialize scalability
  initModel ();

  // apply polarity of BJT
  const char * const type = getPropertyString ("Type");
  pol = !strcmp (type, "pnp") ? -1 : 1;

  // get simulation temperature
  double T = getPropertyDouble ("Temp");

  // initialize starting values
  restartDC ();

  // disable additional base-collector capacitance
  if (deviceEnabled (cbcx)) {
    disableCapacitor (this, cbcx);
  }

  // possibly insert series resistance at emitter
  double Re = getScaledProperty ("Re");
  if (Re != 0.0) {
    // create additional circuit if necessary and reassign nodes
    re = splitResistor (this, re, "Re", "emitter", NODE_E);
    re->setProperty ("R", Re);
    re->setProperty ("Temp", T);
    re->setProperty ("Controlled", getName ());
    re->initDC ();
  }
  // no series resistance at emitter
  else {
    disableResistor (this, re, NODE_E);
  }

  // possibly insert series resistance at collector
  double Rc = getScaledProperty ("Rc");
  if (Rc != 0.0) {
    // create additional circuit if necessary and reassign nodes
    rc = splitResistor (this, rc, "Rc", "collector", NODE_C);
    rc->setProperty ("R", Rc);
    rc->setProperty ("Temp", T);
    rc->setProperty ("Controlled", getName ());
    rc->initDC ();
  }
  // no series resistance at collector
  else {
    disableResistor (this, rc, NODE_C);
  }

  // possibly insert base series resistance
  double Rb  = getScaledProperty ("Rb");
  double Rbm = getScaledProperty ("Rbm");
  if (Rbm <= 0.0) Rbm = Rb; // Rbm defaults to Rb if zero
  if (Rb < Rbm)   Rbm = Rb; // Rbm must be less or equal Rb
  setScaledProperty ("Rbm", Rbm);
  if (Rbm != 0.0) {
    // create additional circuit and reassign nodes
    rb = splitResistor (this, rb, "Rbb", "base", NODE_B);
    rb->setProperty ("R", Rb);
    rb->setProperty ("Temp", T);
    rb->setProperty ("Controlled", getName ());
    rb->initDC ();
  }
  // no series resistance at base
  else {
    disableResistor (this, rb, NODE_B);
    Rbb = 0.0;                 // set this operating point
    setProperty ("Xcjc", 1.0); // other than 1 is senseless here
  }
}

void bjt::restartDC (void) {
  // apply starting values to previous iteration values
  UbePrev = real (getV (NODE_B) - getV (NODE_E)) * pol;
  UbcPrev = real (getV (NODE_B) - getV (NODE_C)) * pol;
}

#define cexState 6 // extra excess phase state

void bjt::calcDC (void) {

  // fetch device model parameters
  double Is   = getScaledProperty ("Is");
  double Nf   = getPropertyDouble ("Nf");
  double Nr   = getPropertyDouble ("Nr");
  double Vaf  = getPropertyDouble ("Vaf");
  double Var  = getPropertyDouble ("Var");
  double Ikf  = getScaledProperty ("Ikf");
  double Ikr  = getScaledProperty ("Ikr");
  double Bf   = getScaledProperty ("Bf");
  double Br   = getScaledProperty ("Br");
  double Ise  = getScaledProperty ("Ise");
  double Isc  = getScaledProperty ("Isc");
  double Ne   = getPropertyDouble ("Ne");
  double Nc   = getPropertyDouble ("Nc");
  double Rb   = getScaledProperty ("Rb");
  double Rbm  = getScaledProperty ("Rbm");
  double Irb  = getScaledProperty ("Irb");
  double T    = getPropertyDouble ("Temp");

  double Ut, Q1, Q2;
  double Iben, Ibcn, Ibei, Ibci, Ibc, gbe, gbc, gtiny;
  double IeqB, IeqC, IeqE, IeqS, UbeCrit, UbcCrit;
  double gm, go;

  // interpret zero as infinity for these model parameters
  Ikf = Ikf > 0 ? 1.0 / Ikf : 0;
  Ikr = Ikr > 0 ? 1.0 / Ikr : 0;
  Vaf = Vaf > 0 ? 1.0 / Vaf : 0;
  Var = Var > 0 ? 1.0 / Var : 0;

  T = celsius2kelvin (T);
  Ut = T * kBoverQ;
  Ube = real (getV (NODE_B) - getV (NODE_E)) * pol;
  Ubc = real (getV (NODE_B) - getV (NODE_C)) * pol;

  // critical voltage necessary for bad start values
  UbeCrit = pnCriticalVoltage (Is, Nf * Ut);
  UbcCrit = pnCriticalVoltage (Is, Nr * Ut);
  UbePrev = Ube = pnVoltage (Ube, UbePrev, Ut * Nf, UbeCrit);
  UbcPrev = Ubc = pnVoltage (Ubc, UbcPrev, Ut * Nr, UbcCrit);

  Uce = Ube - Ubc;

  // base-emitter diodes
  gtiny = Ube < - 10 * Ut * Nf ? (Is + Ise) : 0;
#if 0
  If = pnCurrent (Ube, Is, Ut * Nf);
  Ibei = If / Bf;
  gif = pnConductance (Ube, Is, Ut * Nf);
  gbei = gif / Bf;
  Iben = pnCurrent (Ube, Ise, Ut * Ne);
  gben = pnConductance (Ube, Ise, Ut * Ne);
  Ibe = Ibei + Iben + gtiny * Ube;
  gbe = gbei + gben + gtiny;
#else
  pnJunctionBIP (Ube, Is, Ut * Nf, If, gif);
  Ibei = If / Bf;
  gbei = gif / Bf;
  pnJunctionBIP (Ube, Ise, Ut * Ne, Iben, gben);
  Iben += gtiny * Ube;
  gben += gtiny;
  Ibe = Ibei + Iben;
  gbe = gbei + gben;
#endif

  // base-collector diodes
  gtiny = Ubc < - 10 * Ut * Nr ? (Is + Isc) : 0;
#if 0
  Ir = pnCurrent (Ubc, Is, Ut * Nr);
  Ibci = Ir / Br;
  gir = pnConductance (Ubc, Is, Ut * Nr);
  gbci = gir / Br;
  Ibcn = pnCurrent (Ubc, Isc, Ut * Nc);
  gbcn = pnConductance (Ubc, Isc, Ut * Nc);
  Ibc = Ibci + Ibcn + gtiny * Ubc;
  gbc = gbci + gbcn + gtiny;
#else
  pnJunctionBIP (Ubc, Is, Ut * Nr, Ir, gir);
  Ibci = Ir / Br;
  gbci = gir / Br;
  pnJunctionBIP (Ubc, Isc, Ut * Nc, Ibcn, gbcn);
  Ibcn += gtiny * Ubc;
  gbcn += gtiny;
  Ibc = Ibci + Ibcn;
  gbc = gbci + gbcn;
#endif

  // compute base charge quantities
  Q1 = 1 / (1 - Ubc * Vaf - Ube * Var);
  Q2 = If * Ikf + Ir * Ikr;
  double SArg = 1.0 + 4.0 * Q2;
  double Sqrt = SArg > 0 ? qucs::sqrt (SArg) : 1;
  Qb = Q1 * (1 + Sqrt) / 2;
  dQbdUbe = Q1 * (Qb * Var + gif * Ikf / Sqrt);
  dQbdUbc = Q1 * (Qb * Vaf + gir * Ikr / Sqrt);

  // If and gif will be later used also for the capacitance/charge calculations
  // Values computed from the excess phase routine should be used only
  //   for computing the companion model current and conductance
  double Ifx = If;
  double gifx = gif;
  // during transient analysis only
  if (doTR) {
    // calculate excess phase influence
    Ifx /= Qb;
    excessPhase (cexState, Ifx, gifx);
    Ifx *= Qb;
  }

  // compute transfer current
  It = (Ifx - Ir) / Qb;

  // compute forward and backward transconductance
  gitf = (+gifx - It * dQbdUbe) / Qb;
  gitr = (-gir - It * dQbdUbc) / Qb;

  // compute old SPICE values
  go = -gitr;
  gm = +gitf - go;
  setOperatingPoint ("gm", gm);
  setOperatingPoint ("go", go);

  // calculate current-dependent base resistance
  if (Rbm != 0.0) {
    if (Irb != 0.0) {
      double a, b, z;
      a = (Ibci + Ibcn + Ibei + Iben) / Irb;
      a = std::max (a, NR_TINY); // enforce positive values
      z = (qucs::sqrt (1 + 144 / sqr (pi) * a) - 1) / 24 * sqr (pi) / qucs::sqrt (a);
      b = qucs::tan (z);
      Rbb = Rbm + 3 * (Rb - Rbm) * (b - z) / z / sqr (b);
    }
    else {
      Rbb = Rbm + (Rb - Rbm) / Qb;
    }
    rb->setScaledProperty ("R", Rbb);
    rb->calcDC ();
  }

  // compute autonomic current sources
  IeqB = Ibe - Ube * gbe;
  IeqC = Ibc - Ubc * gbc;
#if NEWSGP
  IeqE = It - Ube * gitf - Ubc * gitr;
#else
  IeqE = It - Ube * gm - Uce * go;
#endif
  IeqS = 0;
  setI (NODE_B, (-IeqB - IeqC) * pol);
  setI (NODE_C, (+IeqC - IeqE - IeqS) * pol);
  setI (NODE_E, (+IeqB + IeqE) * pol);
  setI (NODE_S, (+IeqS) * pol);

  // apply admittance matrix elements
#if NEWSGP
  setY (NODE_B, NODE_B, gbc + gbe);
  setY (NODE_B, NODE_C, -gbc);
  setY (NODE_B, NODE_E, -gbe);
  setY (NODE_B, NODE_S, 0);
  setY (NODE_C, NODE_B, -gbc + gitf + gitr);
  setY (NODE_C, NODE_C, gbc - gitr);
  setY (NODE_C, NODE_E, -gitf);
  setY (NODE_C, NODE_S, 0);
  setY (NODE_E, NODE_B, -gbe - gitf - gitr);
  setY (NODE_E, NODE_C, gitr);
  setY (NODE_E, NODE_E, gbe + gitf);
  setY (NODE_E, NODE_S, 0);
  setY (NODE_S, NODE_B, 0);
  setY (NODE_S, NODE_C, 0);
  setY (NODE_S, NODE_E, 0);
  setY (NODE_S, NODE_S, 0);
#else
  setY (NODE_B, NODE_B, gbc + gbe);
  setY (NODE_B, NODE_C, -gbc);
  setY (NODE_B, NODE_E, -gbe);
  setY (NODE_B, NODE_S, 0);
  setY (NODE_C, NODE_B, -gbc + gm);
  setY (NODE_C, NODE_C, go + gbc);
  setY (NODE_C, NODE_E, -go - gm);
  setY (NODE_C, NODE_S, 0);
  setY (NODE_E, NODE_B, -gbe - gm);
  setY (NODE_E, NODE_C, -go);
  setY (NODE_E, NODE_E, gbe + go + gm);
  setY (NODE_E, NODE_S, 0);
  setY (NODE_S, NODE_B, 0);
  setY (NODE_S, NODE_C, 0);
  setY (NODE_S, NODE_E, 0);
  setY (NODE_S, NODE_S, 0);
#endif
}

void bjt::saveOperatingPoints (void) {
  double Vbe, Vbc;
  Vbe = real (getV (NODE_B) - getV (NODE_E)) * pol;
  Vbc = real (getV (NODE_B) - getV (NODE_C)) * pol;
  Ucs = real (getV (NODE_S) - getV (NODE_C)) * pol;
  setOperatingPoint ("Vbe", Vbe);
  setOperatingPoint ("Vbc", Vbc);
  setOperatingPoint ("Vce", Vbe - Vbc);
  setOperatingPoint ("Vcs", Ucs);
  if (deviceEnabled (cbcx)) {
    Ubx = real (cbcx->getV (NODE_1) - cbcx->getV (NODE_2)) * pol;
    setOperatingPoint ("Vbx", Ubx);
  }
}

void bjt::loadOperatingPoints (void) {
  Ube = getOperatingPoint ("Vbe");
  Ubc = getOperatingPoint ("Vbc");
  Uce = getOperatingPoint ("Vce");
  Ucs = getOperatingPoint ("Vcs");
}

void bjt::calcOperatingPoints (void) {

  // fetch device model parameters
  double Cje0 = getScaledProperty ("Cje");
  double Vje  = getScaledProperty ("Vje");
  double Mje  = getPropertyDouble ("Mje");
  double Cjc0 = getScaledProperty ("Cjc");
  double Vjc  = getScaledProperty ("Vjc");
  double Mjc  = getPropertyDouble ("Mjc");
  double Xcjc = getPropertyDouble ("Xcjc");
  double Cjs0 = getScaledProperty ("Cjs");
  double Vjs  = getScaledProperty ("Vjs");
  double Mjs  = getPropertyDouble ("Mjs");
  double Fc   = getPropertyDouble ("Fc");
  double Vtf  = getPropertyDouble ("Vtf");
  double Tf   = getPropertyDouble ("Tf");
  double Xtf  = getPropertyDouble ("Xtf");
  double Itf  = getScaledProperty ("Itf");
  double Tr   = getPropertyDouble ("Tr");

  double Cbe, Cbci, Cbcx, Ccs;

  // interpret zero as infinity for that model parameter
  Vtf = Vtf > 0 ? 1.0 / Vtf : 0;

  // depletion capacitance of base-emitter diode
  Cbe = pnCapacitance (Ube, Cje0, Vje, Mje, Fc);
  Qbe = pnCharge (Ube, Cje0, Vje, Mje, Fc);

  // diffusion capacitance of base-emitter diode
  if (If != 0.0) {
    double e, Tff, dTffdUbe, dTffdUbc, a;
    a = 1 / (1 + Itf / If);
    e = 2 * qucs::exp (std::min (Ubc * Vtf, 709.0));
    Tff = Tf * (1 + Xtf * sqr (a) * e);
    dTffdUbe = Tf * Xtf * 2 * gif * Itf * cubic (a) / sqr (If) * e;
    Cbe += (If * dTffdUbe + Tff * (gif - If / Qb * dQbdUbe)) / Qb;
    Qbe += If * Tff / Qb;
    dTffdUbc = Tf * Xtf * Vtf * sqr (a) * e;
    dQbedUbc = If / Qb * (dTffdUbc - Tff / Qb * dQbdUbc);
  }

  // depletion and diffusion capacitance of base-collector diode
  Cbci = pnCapacitance (Ubc, Cjc0 * Xcjc, Vjc, Mjc, Fc) + Tr * gir;
  Qbci = pnCharge (Ubc, Cjc0 * Xcjc, Vjc, Mjc, Fc) + Tr * Ir;

  // depletion and diffusion capacitance of external base-collector capacitor
  Cbcx = pnCapacitance (Ubx, Cjc0 * (1 - Xcjc), Vjc, Mjc, Fc);
  Qbcx = pnCharge (Ubx, Cjc0 * (1 - Xcjc), Vjc, Mjc, Fc);

  // depletion capacitance of collector-substrate diode
  Ccs = pnCapacitance (Ucs, Cjs0, Vjs, Mjs);
  Qcs = pnCharge (Ucs, Cjs0, Vjs, Mjs);

  // finally save the operating points
  setOperatingPoint ("Cbe", Cbe);
  setOperatingPoint ("Cbci", Cbci);
  setOperatingPoint ("Cbcx", Cbcx);
  setOperatingPoint ("Ccs", Ccs);
  setOperatingPoint ("gmf", gitf);
  setOperatingPoint ("gmr", gitr);
  setOperatingPoint ("gmu", gbci + gbcn);
  setOperatingPoint ("gpi", gbei + gben);
  setOperatingPoint ("Rbb", Rbb);
  setOperatingPoint ("Ibe", Ibe);
  setOperatingPoint ("Ice", It);
}

void bjt::initSP (void) {
  allocMatrixS ();
  processCbcx ();
  if (deviceEnabled (cbcx)) {
    cbcx->initSP ();
    cbcx->initNoiseSP ();
  }
}

void bjt::processCbcx (void) {
  double Xcjc = getPropertyDouble ("Xcjc");
  double Rbm  = getScaledProperty ("Rbm");
  double Cjc0 = getScaledProperty ("Cjc");

  /* if necessary then insert external capacitance between internal
     collector node and external base node */
  if (Rbm != 0.0 && Cjc0 != 0.0 && Xcjc != 1.0) {
    if (!deviceEnabled (cbcx)) {
      cbcx = splitCapacitor (this, cbcx, "Cbcx", rb->getNode (NODE_1),
			     getNode (NODE_C));
    }
    cbcx->setProperty ("C", getOperatingPoint ("Cbcx"));
  }
  else {
    disableCapacitor (this, cbcx);
  }
}

void bjt::initAC (void) {
  allocMatrixMNA ();
  processCbcx ();
  if (deviceEnabled (cbcx)) {
    cbcx->initAC ();
    cbcx->initNoiseAC ();
  }
}

void bjt::calcAC (double frequency) {
  setMatrixY (calcMatrixY (frequency));
}

void bjt::calcNoiseAC (double frequency) {
  setMatrixN (calcMatrixCy (frequency));
}

#define qbeState 0 // base-emitter charge state
#define cbeState 1 // base-emitter current state
#define qbcState 2 // base-collector charge state
#define cbcState 3 // base-collector current state
#define qcsState 4 // collector-substrate charge state
#define ccsState 5 // collector-substrate current state

#define qbxState 0 // external base-collector charge state
#define cbxState 1 // external base-collector current state

void bjt::initTR (void) {
  setStates (7);
  initDC ();
  doTR = true;

  // handle external base-collector capacitance appropriately
  processCbcx ();
  if (deviceEnabled (cbcx)) {
    cbcx->initTR ();
    cbcx->setProperty ("Controlled", getName ());
  }
}

void bjt::calcTR (double t) {
  calcDC ();
  saveOperatingPoints ();
  loadOperatingPoints ();
  calcOperatingPoints ();

  double Cbe  = getOperatingPoint ("Cbe");
  double Ccs  = getOperatingPoint ("Ccs");
  double Cbci = getOperatingPoint ("Cbci");
  double Cbcx = getOperatingPoint ("Cbcx");

  // handle Rbb and Cbcx appropriately
  if (Rbb != 0.0) {
    rb->setScaledProperty ("R", Rbb);
    rb->calcTR (t);
    if (deviceEnabled (cbcx)) {
      cbcx->clearI ();
      cbcx->clearY ();
      cbcx->transientCapacitance (qbxState, NODE_1, NODE_2, Cbcx, Ubx, Qbcx);
    }
  }

  // usual capacitances
  transientCapacitance (qbeState, NODE_B, NODE_E, Cbe, Ube, Qbe);
  transientCapacitance (qbcState, NODE_B, NODE_C, Cbci, Ubc, Qbci);
  transientCapacitance (qcsState, NODE_S, NODE_C, Ccs, Ucs, Qcs);

  // trans-capacitances
  transientCapacitanceC (NODE_B, NODE_E, NODE_B, NODE_C, dQbedUbc, Ubc);
}

void bjt::excessPhase (int istate, double& i, double& g) {

  // fetch device properties
  double Ptf = getPropertyDouble ("Ptf");
  double Tf = getPropertyDouble ("Tf");
  double td = deg2rad (Ptf) * Tf;

  // return if nothing todo
  if (td == 0.0) return;

  // fill-in current history during initialization
  if (getMode () & MODE_INIT) fillState (istate, i);

  // calculate current coefficients C1, C2 and C3
  double * delta = getDelta ();
  double c3, c2, c1, dn, ra;
  c1 = delta[0] / td;
  c2 = 3 * c1;
  c1 = c2 * c1;
  dn = 1 + c1 + c2;
  c1 = c1 / dn;
  ra = delta[0] / delta[1];
  c2 = (1 + ra + c2) / dn;
  c3 = ra / dn;

  // update and save current, update transconductance
  i = i * c1 + getState (istate, 1) * c2 - getState (istate, 2) * c3;
  setState (istate, i);
  g = g * c1;
}

// properties
PROP_REQ [] = {
  { "Is", PROP_REAL, { 1e-16, PROP_NO_STR }, PROP_POS_RANGE },
  { "Nf", PROP_REAL, { 1, PROP_NO_STR }, PROP_RNGII (0.1, 100) },
  { "Nr", PROP_REAL, { 1, PROP_NO_STR }, PROP_RNGII (0.1, 100) },
  { "Ikf", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Ikr", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Vaf", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Var", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Ise", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Ne", PROP_REAL, { 1.5, PROP_NO_STR }, PROP_RNGII (0.1, 100) },
  { "Isc", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Nc", PROP_REAL, { 2, PROP_NO_STR }, PROP_RNGII (0.1, 100) },
  { "Bf", PROP_REAL, { 100, PROP_NO_STR }, PROP_POS_RANGEX },
  { "Br", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGEX },
  { "Rbm", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Irb", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Cje", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Vje", PROP_REAL, { 0.75, PROP_NO_STR }, PROP_RNGXI (0, 10) },
  { "Mje", PROP_REAL, { 0.33, PROP_NO_STR }, PROP_RNGII (0, 1) },
  { "Cjc", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Vjc", PROP_REAL, { 0.75, PROP_NO_STR }, PROP_RNGXI (0, 10) },
  { "Mjc", PROP_REAL, { 0.33, PROP_NO_STR }, PROP_RNGII (0, 1) },
  { "Xcjc", PROP_REAL, { 1, PROP_NO_STR }, PROP_RNGII (0, 1) },
  { "Cjs", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Vjs", PROP_REAL, { 0.75, PROP_NO_STR }, PROP_RNGXI (0, 10) },
  { "Mjs", PROP_REAL, { 0, PROP_NO_STR }, PROP_RNGII (0, 1) },
  { "Fc", PROP_REAL, { 0.5, PROP_NO_STR }, PROP_RNGII (0, 1) },
  { "Vtf", PROP_REAL, { 0, PROP_NO_STR }, PROP_NO_RANGE },
  { "Tf", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Xtf", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Itf", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Tr", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  PROP_NO_PROP };
PROP_OPT [] = {
  { "Rc", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Re", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Rb", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Kf", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Af", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGE },
  { "Ffe", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGE },
  { "Kb", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
  { "Ab", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGE },
  { "Fb", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGE },
  { "Temp", PROP_REAL, { 26.85, PROP_NO_STR }, PROP_MIN_VAL (K) },
  { "Type", PROP_STR, { PROP_NO_VAL, "npn" }, PROP_RNG_BJT },
  { "Ptf", PROP_REAL, { 0, PROP_NO_STR }, PROP_RNGII (-180, +180) },
  { "Xtb", PROP_REAL, { 0, PROP_NO_STR }, PROP_NO_RANGE },
  { "Xti", PROP_REAL, { 3, PROP_NO_STR }, PROP_POS_RANGE },
  { "Eg", PROP_REAL, { EgSi, PROP_NO_STR }, PROP_POS_RANGE },
  { "Tnom", PROP_REAL, { 26.85, PROP_NO_STR }, PROP_MIN_VAL (K) },
  { "Area", PROP_REAL, { 1, PROP_NO_STR }, PROP_POS_RANGEX },
  PROP_NO_PROP };
struct define_t bjt::cirdef =
  { "BJT", 4, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_NONLINEAR, PROP_DEF };
