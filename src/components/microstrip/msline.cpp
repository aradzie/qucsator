/*
 * msline.cpp - microstrip transmission line class implementation
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
#include "substrate.h"
#include "msline.h"

using namespace qucs;

msline::msline () : circuit (2) {
  alpha = beta = zl = ereff = 0;
  type = CIR_MSLINE;
}

void msline::calcNoiseSP (double) {
  double l = getPropertyDouble ("L");
  if (l < 0) return;
  // calculate noise using Bosma's theorem
  double T = getPropertyDouble ("Temp");
  matrix s = getMatrixS ();
  matrix e = eye (getSize ());
  setMatrixN (celsius2kelvin (T) / T0 * (e - s * transpose (conj (s))));
}

void msline::calcPropagation (double frequency) {

  /* how to get properties of this component, e.g. L, W */
  double W = getPropertyDouble ("W");
  const char * SModel = getPropertyString ("Model");
  const char * DModel = getPropertyString ("DispModel");

  /* how to get properties of the substrate, e.g. Er, H */
  substrate * subst = getSubstrate ();
  double er    = subst->getPropertyDouble ("er");
  double h     = subst->getPropertyDouble ("h");
  double t     = subst->getPropertyDouble ("t");
  double tand  = subst->getPropertyDouble ("tand");
  double rho   = subst->getPropertyDouble ("rho");
  double D     = subst->getPropertyDouble ("D");

  /* local variables */
  double ac, ad;
  double ZlEff, ErEff, WEff, ZlEffFreq, ErEffFreq;

  // quasi-static effective dielectric constant of substrate + line and
  // the impedance of the microstrip line
  analyseQuasiStatic (W, h, t, er, SModel, ZlEff, ErEff, WEff);

  // analyse dispersion of Zl and Er (use WEff here?)
  analyseDispersion (W, h, er, ZlEff, ErEff, frequency, DModel,
		     ZlEffFreq, ErEffFreq);

  // analyse losses of line
  analyseLoss (W, t, er, rho, D, tand, ZlEff, ZlEff, ErEff,
	       frequency, "Hammerstad", ac, ad);

  // calculate propagation constants and reference impedance
  zl    = ZlEffFreq;
  ereff = ErEffFreq;
  alpha = ac + ad;
  beta  = qucs::sqrt (ErEffFreq) * 2 * pi * frequency / C0;
}

void msline::calcSP (double frequency) {
  double l = getPropertyDouble ("L");

  // calculate propagation constants
  calcPropagation (frequency);

  // calculate S-parameters
  double z = zl / z0;
  double y = 1 / z;
  nr_complex_t g = nr_complex_t (alpha, beta);
  nr_complex_t n = 2.0 * cosh (g * l) + (z + y) * qucs::sinh (g * l);
  nr_complex_t s11 = (z - y) * qucs::sinh (g * l) / n;
  nr_complex_t s21 = 2.0 / n;
  setS (NODE_1, NODE_1, s11); setS (NODE_2, NODE_2, s11);
  setS (NODE_1, NODE_2, s21); setS (NODE_2, NODE_1, s21);
}

void msline::saveCharacteristics (double) {
  setCharacteristic ("Zl", zl);
  setCharacteristic ("Er", ereff);
}

/* This function calculates the quasi-static impedance of a microstrip
   line, the value of the effective dielectric constant and the
   effective width due to the finite conductor thickness for the given
   microstrip line and substrate properties. */
void msline::analyseQuasiStatic (double W, double h, double t,
				 double er, const char * const Model,
				 double& ZlEff, double& ErEff,
				 double& WEff) {

  double z, e;

  // default values
  e = er;
  z = z0;
  WEff = W;

  // WHEELER
  if (!strcmp (Model, "Wheeler")) {
    double a, b, c, d, x, dW1, dWr, Wr;

    // compute strip thickness effect
    if (t != 0) {
      dW1 = t / pi * qucs::log (4 * euler / qucs::sqrt (sqr (t / h) +
					    sqr (one_over_pi / (W / t + 1.10))));
    }
    else dW1 = 0;
    dWr = (1 + 1 / er) / 2 * dW1;
    Wr  = WEff = W + dWr;

    // compute characteristic impedance
    if (W / h < 3.3) {
      c = qucs::log (4 * h / Wr + qucs::sqrt (sqr (4 * h / Wr) + 2));
      b = (er - 1) / (er + 1) / 2 * (qucs::log (pi_over_2) + qucs::log (2 * two_over_pi) / er);
      z = (c - b) * Z0 / pi / qucs::sqrt (2 * (er + 1));
    }
    else {
      c = 1 + qucs::log (pi_over_2) + qucs::log (Wr / h / 2 + 0.94);
      d = one_over_pi / 2 * (1 + qucs::log (sqr (pi) / 16)) * (er - 1) / sqr (er);
      x = 2 * ln2 / pi + Wr / h / 2 + (er + 1) / 2 / pi / er * c + d;
      z = Z0 / 2 / x / qucs::sqrt (er);
    }

    // compute effective dielectric constant
    if (W / h < 1.3) {
      a = qucs::log (8 * h / Wr) + sqr (Wr / h) / 32;
      b = (er - 1) / (er + 1) / 2 * (qucs::log (pi_over_2) + qucs::log (2 * two_over_pi) / er);
      e = (er + 1) / 2 * sqr (a / (a - b));
    }
    else {
      a = (er - 1) / 2 / pi / er * (qucs::log (2.1349 * Wr / h + 4.0137) -
				      0.5169 / er);
      b = Wr / h / 2 + one_over_pi * qucs::log (8.5397 * Wr / h + 16.0547);
      e = er * sqr ((b - a) / b);
    }
  }
  // SCHNEIDER
  else if (!strcmp (Model, "Schneider")) {

    double dW = 0, u = W / h;

    // consider strip thickness equations
    if (t != 0 && t < W / 2) {
      double arg = (u < one_over_pi / 2) ? 2 * pi * W / t : h / t;
      dW = t / pi * (1 + qucs::log (2 * arg));
      if (t / dW >= 0.75) dW = 0;
    }
    WEff = W + dW; u = WEff / h;

    // effective dielectric constant
    e = (er + 1) / 2 + (er - 1) / 2 / qucs::sqrt (1 + 10 / u);

    // characteristic impedance
    if (u < 1.0) {
      z = one_over_pi / 2 * qucs::log (8 / u + u / 4);
    }
    else {
      z = 1 / (u + 2.42 - 0.44 / u + qucs::pow ((1. - 1. / u), 6.));
    }
    z = Z0 * z / qucs::sqrt (e);
  }
  // HAMMERSTAD and JENSEN
  else if (!strcmp (Model, "Hammerstad")) {
    double a, b, du1, du, u, ur, u1, zr, z1;

    u = W / h; // normalized width
    t = t / h; // normalized thickness

    // compute strip thickness effect
    if (t != 0) {
      du1 = t / pi * qucs::log (1 + 4 * euler / t / sqr (coth (qucs::sqrt (6.517 * u))));
    }
    else du1 = 0;
    du = du1 * (1 + sech (qucs::sqrt (er - 1))) / 2;
    u1 = u + du1;
    ur = u + du;
    WEff = ur * h;

    // compute impedances for homogeneous medium
    Hammerstad_zl (ur, zr);
    Hammerstad_zl (u1, z1);

    // compute effective dielectric constant
    Hammerstad_ab (ur, er, a, b);
    Hammerstad_er (ur, er, a, b, e);

    // compute final characteristic impedance and dielectric constant
    // including strip thickness effects
    z = zr / qucs::sqrt (e);
    e = e * sqr (z1 / zr);
  }

  ZlEff = z;
  ErEff = e;
}

/* This function calculates the frequency dependent value of the
   effective dielectric constant and the microstrip line impedance for
   the given frequency. */
void msline::analyseDispersion (double W, double h, double er,
				double ZlEff, double ErEff,
				double frequency, const char * const Model,
				double& ZlEffFreq,
				double& ErEffFreq) {

  double e, z;

  // default values
  z = ZlEffFreq = ZlEff;
  e = ErEffFreq = ErEff;

  // GETSINGER
  if (!strcmp (Model, "Getsinger")) {
    Getsinger_disp (h, er, ErEff, ZlEff, frequency, e, z);
  }
  // SCHNEIDER
  else if (!strcmp (Model, "Schneider")) {
    double k, f;
    k = qucs::sqrt (ErEff / er);
    f = 4 * h * frequency / C0 * qucs::sqrt (er - 1);
    f = sqr (f);
    e = ErEff * sqr ((1 + f) / (1 + k * f));
    z = ZlEff * qucs::sqrt (ErEff / e);
  }
  // YAMASHITA
  else if (!strcmp (Model, "Yamashita")) {
    double k, f;
    k = qucs::sqrt (er / ErEff);
    f = 4 * h * frequency / C0 * qucs::sqrt (er - 1) *
      (0.5 + sqr (1 + 2 * qucs::log10 (1 + W / h)));
    e = ErEff * sqr ((1 + k * qucs::pow (f, 1.5) / 4) / (1 + qucs::pow (f, 1.5) / 4));
  }
  // KOBAYASHI
  else if (!strcmp (Model, "Kobayashi")) {
    double n, no, nc, fh, fk;
    fk = C0 * qucs::atan (er * qucs::sqrt ((ErEff - 1) / (er - ErEff))) /
      (2 * pi * h * qucs::sqrt (er - ErEff));
    fh = fk / (0.75 + (0.75 - 0.332 / qucs::pow (er, 1.73)) * W / h);
    no = 1 + 1 / (1 + qucs::sqrt (W / h)) + 0.32 * cubic (1 / (1 + qucs::sqrt (W / h)));
    if (W / h < 0.7) {
      nc = 1 + 1.4 / (1 + W / h) * (0.15 - 0.235 *
				    qucs::exp (-0.45 * frequency / fh));
    }
    else nc = 1;
    n = no * nc < 2.32 ? no * nc : 2.32;
    e = er - (er - ErEff) / (1 + qucs::pow (frequency / fh, n));
  }
  // PRAMANICK and BHARTIA
  else if (!strcmp (Model, "Pramanick")) {
    double Weff, We, f;
    f = 2 * MU0 * h * frequency * qucs::sqrt (ErEff / er) / ZlEff;
    e = er - (er - ErEff) / (1 + sqr (f));
    Weff = Z0 * h / ZlEff / qucs::sqrt (ErEff);
    We = W + (Weff - W) / (1 + sqr (f));
    z = Z0 * h / We / qucs::sqrt (e);
  }
  // HAMMERSTAD and JENSEN
  else if (!strcmp (Model, "Hammerstad")) {
    double f, g;
    g = sqr (pi) / 12 * (er - 1) / ErEff * qucs::sqrt (2 * pi * ZlEff / Z0);
    f = 2 * MU0 * h * frequency / ZlEff;
    e = er - (er - ErEff) / (1 + g * sqr (f));
    z = ZlEff * qucs::sqrt (ErEff / e) * (e - 1) / (ErEff - 1);
  }
  // KIRSCHNING and JANSEN
  else if (!strcmp (Model, "Kirschning")) {
    double r17, u  = W / h, fn = frequency * h / 1e6;

    // dispersion of dielectric constant
    Kirschning_er (u, fn, er, ErEff, e);

    // dispersion of characteristic impedance
    Kirschning_zl (u, fn, er, ErEff, e, ZlEff, r17, z);
  }

  ZlEffFreq = z;
  ErEffFreq = e;
}

/* Computes the exponent factors a(u) and b(er) used within the
   effective relative dielectric constant calculations for single and
   coupled microstrip lines by Hammerstad and Jensen. */
void msline::Hammerstad_ab (double u, double er, double& a,
			    double& b) {
  a = 1 + qucs::log ((quadr (u) + sqr (u / 52)) / (quadr (u) + 0.432)) / 49 +
    qucs::log (1 + cubic (u / 18.1)) / 18.7;
  b = 0.564 * qucs::pow ((er - 0.9) / (er + 3), 0.053);
}

/* The function computes the effective dielectric constant of a single
   microstrip.  The equation is used in single and coupled microstrip
   calculations. */
void msline::Hammerstad_er (double u, double er, double a,
			    double b, double& e) {
  e = (er + 1) / 2 + (er - 1) / 2 * qucs::pow (1 + 10 / u, -a * b);
}

/* This function computes the characteristic impedance of single
   microstrip line based upon the given width-height ratio.  The
   equation is used in single and coupled microstrip calculations as
   well. */
void msline::Hammerstad_zl (double u, double& zl) {
  double fu = 6 + (2 * pi - 6) * qucs::exp (- qucs::pow (30.666 / u, 0.7528));
  zl = Z0 / 2 / pi * qucs::log (fu / u + qucs::sqrt (1 + sqr (2 / u)));
}

/* Calculates dispersion effects for effective dielectric constant and
   characteristic impedance as defined by Getsinger (for single and
   coupled microstrips). */
void msline::Getsinger_disp (double h, double er, double ErEff,
			     double ZlEff, double frequency,
			     double& e, double& z) {
  double g, f, d;
  g = 0.6 + 0.009 * ZlEff;
  f = frequency * 2 * MU0 * h / ZlEff;
  e = er - (er - ErEff) / (1 + g * sqr (f));
  d = (er - e) * (e - ErEff) / e / (er - ErEff);
  z = ZlEff * qucs::sqrt (e / ErEff) / (1 + d);  // group delay model
}

/* This function computes the dispersion of the effective dielectric
   constant of a single microstrip line.  It is defined in a separate
   function because it is used within the coupled microstrip lines as
   well. */
void msline::Kirschning_er (double u, double fn, double er,
			    double ErEff, double& ErEffFreq) {
  double p, p1, p2, p3, p4;
  p1 = 0.27488 + (0.6315 + 0.525 / qucs::pow (1. + 0.0157 * fn, 20.)) * u -
    0.065683 * qucs::exp (-8.7513 * u);
  p2 = 0.33622 * (1 - qucs::exp (-0.03442 * er));
  p3 = 0.0363 * qucs::exp (-4.6 * u) * (1 - qucs::exp (- qucs::pow (fn / 38.7, 4.97)));
  p4 = 1 + 2.751 * (1 - qucs::exp (- qucs::pow (er / 15.916, 8.)));
  p  = p1 * p2 * qucs::pow ((0.1844 + p3 * p4) * fn, 1.5763);
  ErEffFreq  = er - (er - ErEff) / (1 + p);
}

/* Computes dispersion effects of characteristic impedance of a single
   microstrip line according to Kirschning and Jansen.  Also used in
   coupled microstrip lines calculations. */
void msline::Kirschning_zl (double u, double fn, double er,
			    double ErEff, double ErEffFreq,
			    double ZlEff, double& r17,
			    double& ZlEffFreq) {
  double r1, r2, r3, r4, r5, r6, r7, r8, r9, r10;
  double r11, r12, r13, r14, r15, r16;
  r1 = 0.03891 * qucs::pow (er, 1.4);
  r2 = 0.267 * qucs::pow (u, 7.);
  r3 = 4.766 * qucs::exp (-3.228 * qucs::pow (u, 0.641));
  r4 = 0.016 + qucs::pow (0.0514 * er, 4.524);
  r5 = qucs::pow (fn / 28.843, 12.);
  r6 = 22.20 * qucs::pow (u, 1.92);
  r7 = 1.206 - 0.3144 * qucs::exp (-r1) * (1 - qucs::exp (-r2));
  r8 = 1 + 1.275 * (1 - qucs::exp (-0.004625 * r3 *
			     qucs::pow (er, 1.674) * qucs::pow (fn / 18.365, 2.745)));
  r9 = 5.086 * r4 * r5 / (0.3838 + 0.386 * r4) *
    qucs::exp (-r6) / (1 + 1.2992 * r5) *
    qucs::pow (er - 1., 6.) / (1 + 10 * qucs::pow (er - 1., 6.));
  r10 = 0.00044 * qucs::pow (er, 2.136) + 0.0184;
  r11 = qucs::pow (fn / 19.47, 6.) / (1 + 0.0962 * qucs::pow (fn / 19.47, 6.));
  r12 = 1 / (1 + 0.00245 * sqr (u));
  r13 = 0.9408 * qucs::pow (ErEffFreq, r8) - 0.9603;
  r14 = (0.9408 - r9) * qucs::pow (ErEff, r8) - 0.9603;
  r15 = 0.707 * r10 * qucs::pow (fn / 12.3, 1.097);
  r16 = 1 + 0.0503 * sqr (er) * r11 * (1 - qucs::exp (- qucs::pow (u / 15., 6.)));
  r17 = r7 * (1 - 1.1241 * r12 / r16 *
	      qucs::exp (-0.026 * qucs::pow (fn, 1.15656) - r15));
  ZlEffFreq = ZlEff * qucs::pow (r13 / r14, r17);
}

/* The function calculates the conductor and dielectric losses of a
   single microstrip line. */
void msline::analyseLoss (double W, double t, double er,
			  double rho, double D, double tand,
			  double ZlEff1, double ZlEff2,
			  double ErEff,
			  double frequency, const char * Model,
			  double& ac, double& ad) {
  ac = ad = 0;

  // HAMMERSTAD and JENSEN
  if (!strcmp (Model, "Hammerstad")) {
    double Rs, ds, l0, Kr, Ki;

    // conductor losses
    if (t != 0.0) {
      Rs = qucs::sqrt (pi * frequency * MU0 * rho); // skin resistance
      ds = rho / Rs;                            // skin depth
      // valid for t > 3 * ds
      if (t < 3 * ds) {
	logprint (LOG_ERROR,
		  "WARNING: conductor loss calculation invalid for line "
		  "height t (%g) < 3 * skin depth (%g)\n", t, 3 * ds);
      }
      // current distribution factor
      Ki = qucs::exp (-1.2 * qucs::pow ((ZlEff1 + ZlEff2) / 2 / Z0, 0.7));
      // D is RMS surface roughness
      Kr = 1 + two_over_pi * qucs::atan (1.4 * sqr (D / ds));
      ac = Rs / (ZlEff1 * W) * Ki * Kr;
    }

    // dielectric losses
    l0 = C0 / frequency;
    ad = pi * er / (er - 1) * (ErEff - 1) / qucs::sqrt (ErEff) * tand / l0;
  }
}

void msline::initDC (void) {
  double l     = getPropertyDouble ("L");
  double W     = getPropertyDouble ("W");
  substrate * subst = getSubstrate ();
  double t     = subst->getPropertyDouble ("t");
  double rho   = subst->getPropertyDouble ("rho");

  if (t != 0.0 && rho != 0.0 && l != 0.0) {
    // tiny resistance
    double g = t * W / rho / l;
    setVoltageSources (0);
    allocMatrixMNA ();
    setY (NODE_1, NODE_1, +g); setY (NODE_2, NODE_2, +g);
    setY (NODE_1, NODE_2, -g); setY (NODE_2, NODE_1, -g);
  }
  else {
    // a DC short (voltage source V = 0 volts)
    setVoltageSources (1);
    setInternalVoltageSource (1);
    allocMatrixMNA ();
    clearY ();
    voltageSource (VSRC_1, NODE_1, NODE_2);
  }
}

void msline::initAC (void) {
  setVoltageSources (0);
  allocMatrixMNA ();
}

void msline::calcAC (double frequency) {
  double l = getPropertyDouble ("L");

  // calculate propagation constants
  calcPropagation (frequency);

  // calculate Y-parameters
  nr_complex_t g = nr_complex_t (alpha, beta);
  nr_complex_t y11 = coth (g * l) / zl;
  nr_complex_t y21 = -cosech (g * l) / zl;
  setY (NODE_1, NODE_1, y11); setY (NODE_2, NODE_2, y11);
  setY (NODE_1, NODE_2, y21); setY (NODE_2, NODE_1, y21);
}

void msline::calcNoiseAC (double) {
  double l = getPropertyDouble ("L");
  if (l < 0) return;
  // calculate noise using Bosma's theorem
  double T = getPropertyDouble ("Temp");
  setMatrixN (4 * celsius2kelvin (T) / T0 * real (getMatrixY ()));
}

// properties
PROP_REQ [] = {
  { "W", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "L", PROP_REAL, { 10e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "Subst", PROP_STR, { PROP_NO_VAL, "Subst1" }, PROP_NO_RANGE },
  { "DispModel", PROP_STR, { PROP_NO_VAL, "Kirschning" }, PROP_RNG_DIS },
  { "Model", PROP_STR, { PROP_NO_VAL, "Hammerstad" }, PROP_RNG_MOD },
  PROP_NO_PROP };
PROP_OPT [] = {
  { "Temp", PROP_REAL, { 26.85, PROP_NO_STR }, PROP_MIN_VAL (K) },
  PROP_NO_PROP };
struct define_t msline::cirdef =
  { "MLIN", 2, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
