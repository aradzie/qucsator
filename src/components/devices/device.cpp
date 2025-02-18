/*
 * device.cpp - device class implementation
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
 */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

#include "complex.h"
#include "object.h"
#include "node.h"
#include "circuit.h"
#include "net.h"
#include "constants.h"
#include "device.h"
#include "netdefs.h"
#include "resistor.h"
#include "capacitor.h"

using namespace qucs;
using namespace qucs::device;

/* This function can be used to create an extra resistor circuit.  If
   the 'res' argument is NULL then the new circuit is created, the
   nodes get re-arranged and it is inserted into the given
   netlist. The given arguments can be explained as follows.
   base:     calling circuit (this)
   res:      additional resistor circuit (can be NULL)
   c:        name of the additional circuit
   n:        name of the inserted (internal) node
   internal: number of new internal node (the original external node) */
circuit * device::splitResistor (circuit * base, circuit * res,
				 const char * c, const char * n,
				 int internal) {
  if (res == NULL) {
    res = new resistor ();
    const std::string &name = circuit::createInternal (c, base->getName ());
    const std::string &node = circuit::createInternal (n, base->getName ());
    res->setName (name);
    res->setNode (0, base->getNode(internal)->getName ());
    res->setNode (1, node, 1);
    base->getNet()->insertCircuit (res);
  }
  base->setNode (internal, res->getNode(1)->getName (), 1);
  return res;
}

/* This function is the counterpart of the above routine.  It removes
   the resistor circuit from the netlist and re-assigns the original
   node. */
void device::disableResistor (circuit * base, circuit * res, int internal) {
  if (res != NULL) {
    base->getNet()->removeCircuit (res, 0);
    base->setNode (internal, res->getNode(1)->getName (), 0);
  }
}

/* This function creates a new capacitor circuit if the given one is
   not NULL.  The new circuit is connected between the given nodes and
   a name is applied based upon the parents (base) name and the given
   name 'c'.  The circuit is then put into the netlist. */
circuit * device::splitCapacitor (circuit * base, circuit * cap,
				  const char * c, node * n1, node * n2) {
  if (cap == NULL) {
    cap = new capacitor ();
    const std::string &name = circuit::createInternal (c, base->getName ());
    cap->setName (name);
    cap->setNode (0, n1->getName ());
    cap->setNode (1, n2->getName ());
  }
  base->getNet()->insertCircuit (cap);
  return cap;
}

// The function removes the given capacitor circuit from the netlist.
void device::disableCapacitor (circuit * base, circuit * cap) {
  if (cap != NULL) {
    base->getNet()->removeCircuit (cap, 0);
  }
}

/* This function checks whether the given circuit object exists and is
   chained within the current netlist.  It returns non-zero if so and
   zero otherwise. */
int device::deviceEnabled (circuit * c) {
  if (c != NULL && c->isEnabled ())
    return 1;
  return 0;
}

/* The function limits the forward pn-voltage for each DC iteration in
   order to avoid numerical overflows and thereby improve the
   convergence. */
double device::pnVoltage (double Ud, double Uold,
			       double Ut, double Ucrit) {
  double arg;
  if (Ud > Ucrit && fabs (Ud - Uold) > 2 * Ut) {
    if (Uold > 0) {
      arg = (Ud - Uold) / Ut;
      if (arg > 0)
	Ud = Uold + Ut * (2 + log (arg - 2));
      else
	Ud = Uold - Ut * (2 + log (2 - arg));
    }
    else Ud = Uold < 0 ? Ut * log (Ud / Ut) : Ucrit;
  }
  else {
    if (Ud < 0) {
      arg = Uold > 0 ? -1 - Uold : 2 * Uold - 1;
      if (Ud < arg) Ud = arg;
    }
  }
  return Ud;
}

// Computes current and its derivative for a MOS pn-junction.
void device::pnJunctionMOS (double Upn, double Iss, double Ute,
			    double& I, double& g) {
  if (Upn <= 0) {
    g = Iss / Ute;
    I = g * Upn;
  }
  else {
    double e = exp (std::min (Upn / Ute, 709.0));
    I = Iss * (e - 1);
    g = Iss * e / Ute;
  }
}

// Computes current and its derivative for a bipolar pn-junction.
void device::pnJunctionBIP (double Upn, double Iss, double Ute,
			    double& I, double& g) {
  if (Upn < -3 * Ute) {
    double a = 3 * Ute / (Upn * euler);
    a = cubic (a);
    I = -Iss * (1 + a);
    g = +Iss * 3 * a / Upn;
  }
  else {
    double e = exp (std::min (Upn / Ute, 709.0));
    I = Iss * (e - 1);
    g = Iss * e / Ute;
  }
}

// The function computes the exponential pn-junction current.
double
device::pnCurrent (double Upn, double Iss, double Ute) {
  return Iss * (exp (std::min (Upn / Ute, 709.0)) - 1);
}

// The function computes the exponential pn-junction current's derivative.
double
device::pnConductance (double Upn, double Iss, double Ute) {
  return Iss * exp (std::min (Upn / Ute, 709.0)) / Ute;
}

// Computes pn-junction depletion capacitance.
double
device::pnCapacitance (double Uj, double Cj, double Vj,
		       double Mj, double Fc) {
  double c;
  if (Uj <= Fc * Vj)
    c = Cj * exp (-Mj * log (1 - Uj / Vj));
  else
    c = Cj * exp (-Mj * log (1 - Fc)) *
      (1 + Mj * (Uj - Fc * Vj) / Vj / (1 - Fc));
  return c;
}

// Computes pn-junction depletion charge.
double device::pnCharge (double Uj, double Cj, double Vj,
			      double Mj, double Fc) {
  double q, a, b;
  if (Uj <= Fc * Vj) {
    a = 1 - Uj / Vj;
    b = exp ((1 - Mj) * log (a));
    q = Cj * Vj / (1 - Mj) * (1 - b);
  }
  else {
#if 0
    a = 1 - Fc;
    b = exp ((1 - Mj) * log (a));
    a = exp ((1 + Mj) * log (a));
    double c = 1 - Fc * (1 + Mj);
    double d = Fc * Vj;
    double e = Vj * (1 - b) / (1 - Mj);
    q = Cj * (e + (c * (Uj - d) + Mj / 2 / Vj * (sqr (Uj) - sqr (d))) / a);
#else /* this variant is numerically more stable */
    a = 1 - Fc;
    b = exp (-Mj * log (a));
    double f = Fc * Vj;
    double c = Cj * (1 - Fc * (1 + Mj)) * b / a;
    double d = Cj * Mj * b / a / Vj;
    double e = Cj * Vj * (1 - a * b) / (1 - Mj) - d / 2 * f * f - f * c;
    q = e + Uj * (c + Uj * d / 2);
#endif
  }
  return q;
}

/* This function computes the pn-junction depletion capacitance with
   no linearization factor given. */
double
device::pnCapacitance (double Uj, double Cj, double Vj,
		       double Mj) {
  double c;
  if (Uj <= 0)
    c = Cj * exp (-Mj * log (1 - Uj / Vj));
  else
    c = Cj * (1 + Mj * Uj / Vj);
  return c;
}

/* This function computes the pn-junction depletion charge with no
   linearization factor given. */
double device::pnCharge (double Uj, double Cj, double Vj,
			      double Mj) {
  double q;
  if (Uj <= 0)
    q = Cj * Vj / (1 - Mj) * (1 - exp ((1 - Mj) * log (1 - Uj / Vj)));
  else
    q = Cj * Uj * (1 + Mj * Uj / 2 / Vj);
  return q;
}

// Compute critical voltage of pn-junction.
double device::pnCriticalVoltage (double Iss, double Ute) {
  return Ute * log (Ute / sqrt2 / Iss);
}

/* The function limits the forward fet-voltage for each DC iteration
   in order to avoid numerical overflows and thereby improve the
   convergence. */
double
device::fetVoltage (double Ufet, double Uold, double Uth) {
  double Utsthi = fabs (2 * (Uold - Uth)) + 2.0;
  double Utstlo = Utsthi / 2;
  double Utox   = Uth + 3.5;
  double DeltaU = Ufet - Uold;

  if (Uold >= Uth) { /* FET is on */
    if (Uold >= Utox) {
      if (DeltaU <= 0) { /* going off */
	if (Ufet >= Utox) {
	  if (-DeltaU > Utstlo) {
	    Ufet = Uold - Utstlo;
	  }
	} else {
	  Ufet = std::max (Ufet, Uth + 2);
	}
      } else { /* staying on */
	if (DeltaU >= Utsthi) {
	  Ufet = Uold + Utsthi;
	}
      }
    } else { /* middle region */
      if (DeltaU <= 0) { /* decreasing */
	Ufet = std::max (Ufet, Uth - 0.5);
      } else { /* increasing */
	Ufet = std::min (Ufet, Uth + 4);
      }
    }
  } else { /* FET is off */
    if (DeltaU <= 0) { /* staying off */
      if (-DeltaU > Utsthi) {
	Ufet = Uold - Utsthi;
      }
    } else { /* going on */
      if (Ufet <= Uth + 0.5) {
	if (DeltaU > Utstlo) {
	  Ufet = Uold + Utstlo;
	}
      } else {
	Ufet = Uth + 0.5;
      }
    }
  }
  return Ufet;
}

/* The function limits the drain-source voltage for each DC iteration
   in order to avoid numerical overflows and thereby improve the
   convergence. */
double device::fetVoltageDS (double Ufet, double Uold) {
  if (Uold >= 3.5) {
    if (Ufet > Uold) {
      Ufet = std::min (Ufet, 3 * Uold + 2);
    } else if (Ufet < 3.5) {
      Ufet = std::max (Ufet, 2.0);
    }
  } else {
    if (Ufet > Uold) {
      Ufet = std::min (Ufet, 4.0);
    } else {
      Ufet = std::max (Ufet, -0.5);
    }
  }
  return Ufet;
}

/* This function calculates the overlap capacitance for MOS based upon
   the given voltages, surface potential and the zero-bias oxide
   capacitance. */
void
device::fetCapacitanceMeyer (double Ugs, double Ugd, double Uth,
			     double Udsat, double Phi,
			     double Cox, double& Cgs,
			     double& Cgd, double& Cgb) {

  double Utst = Ugs - Uth;
  if (Utst <= -Phi) { // cutoff regions
    Cgb = Cox;
    Cgs = 0;
    Cgd = 0;
  } else if (Utst <= -Phi / 2) {
    Cgb = -Utst * Cox / Phi;
    Cgs = 0;
    Cgd = 0;
  } else if (Utst <= 0) { // depletion region
    Cgb = -Utst * Cox / Phi;
    Cgs = Utst * Cox * 4 / 3 / Phi + 2 * Cox / 3;
    Cgd = 0;
  } else {
    Cgb = 0;
    double Uds = Ugs - Ugd;
    if (Udsat <= Uds) { // saturation region
      Cgs = 2 * Cox / 3;
      Cgd = 0;
    } else { // linear region
      double Sqr1 = sqr (Udsat - Uds);
      double Sqr2 = sqr (2 * Udsat - Uds);
      Cgs = Cox * (1 - Sqr1 / Sqr2) * 2 / 3;
      Cgd = Cox * (1 - Udsat * Udsat / Sqr2) * 2 / 3;
    }
  }
}

// Computes temperature dependency of energy bandgap.
double device::Egap (double T, double Eg0) {
  double a = 7.02e-4;
  double b = 1108;
  return Eg0 - (a * sqr (T)) / (T + b);
}

// Computes temperature dependency of intrinsic density.
double device::intrinsicDensity (double T, double Eg0) {
  double TR = 300.00;
  double E1 = Egap (TR, Eg0);
  double E2 = Egap (T, Eg0);
  double NI = NiSi / 1e6;
  return  NI * exp (1.5 * log (T / TR) + (E1 / TR - E2 / T) / kBoverQ / 2);
}

// Calculates temperature dependence for saturation current.
double
device::pnCurrent_T (double T1, double T2, double Is,
		     double Eg, double N,  double Xti) {
  double Vt, TR;
  TR = T2 / T1;
  Vt = T2 * kBoverQ;
  return Is * exp (Xti / N * log (TR) - Eg / N / Vt * (1 - TR));
}

// Calculates temperature dependence for junction potential.
double
device::pnPotential_T (double T1, double T2, double Vj,
		       double Eg0) {
  double Vt, TR, E1, E2;
  TR = T2 / T1;
  E1 = Egap (T1, Eg0);
  E2 = Egap (T2, Eg0);
  Vt = T2 * kBoverQ;
  return TR * Vj - 3 * Vt * log (TR) - (TR * E1 - E2);
}

// Calculates temperature dependence for junction capacitance.
double
device::pnCapacitance_T (double T1, double T2, double M,
			 double VR, double Cj) {
  return Cj * pnCapacitance_F (T1, T2, M, VR);
}

// Calculates temperature dependence for junction capacitance.
double
device::pnCapacitance_F (double T1, double T2, double M,
			 double VR) {
  double DT = T2 - T1;
  return 1 + M * (4e-4 * DT - VR + 1);
}
