/*
 * transient.cpp - transient helper class implementation
 *
 * Copyright (C) 2004, 2006 Stefan Jahn <stefan@lkcc.org>
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

#include "transient.h"
#include "eqnsys.h"
#include "integrator.h"
#include "tmatrix.h"
#include "tvector.h"

namespace qucs {

constexpr bool FIXED_COEFF = false;

using namespace transient;

/* Computes the integration coefficient for numerical integration methods.
 * Supported methods are: Gear (order 1-6), Trapezoidal, backward Euler
 * and Adams-Moulton (order 1-6). */
void transient::calcCorrectorCoeff(const int method, const int order, const double *deltas,
                                   double *coefficients) {
  tmatrix<double> A(order + 1);
  tvector<double> x(order + 1);
  tvector<double> b(order + 1);
  eqnsys<double> e;
  e.setAlgo(ALGO_LU_DECOMPOSITION);

  switch (method) {
  case INTEGRATOR_EULER: {
    // BACKWARD EULER
    coefficients[0] = 1 / deltas[0];
    coefficients[1] = -1 / deltas[0];
    break;
  }
  case INTEGRATOR_TRAPEZOIDAL: {
    // TRAPEZOIDAL (bilinear)
    coefficients[0] = 2 / deltas[0];
    coefficients[1] = -2 / deltas[0];
    break;
  }
  case INTEGRATOR_GEAR: {
    // GEAR order 1 to 6
    if constexpr (FIXED_COEFF) {
      // right hand side vector
      for (int i = 0; i < order + 1; i++) {
        b.set(i, 1);
      }
      for (int i = 1; i < order + 1; i++) {
        A.set(i, 0, i); // first column
        A.set(0, i, 1); // first row
      }
      for (int c = 1; c <= order - 1; c++) {
        double entry = -c;
        for (int r = 1; r <= order; r++) {
          A.set(r, c + 1, entry);
          entry *= -c;
        }
      }
      e.passEquationSys(&A, &x, &b);
      e.solve();
      // vector x consists of b_{-1}, a_{0}, a_{1} ... a_{k-1} right here
      double k = x.get(0);
      coefficients[0] = 1 / deltas[0] / k;
      for (int i = 1; i <= order; i++) {
        coefficients[i] = -1 / deltas[0] / k * x.get(i);
      }
    } else {
      // right hand side vector
      b.set(1, -1 / deltas[0]);
      // first row
      for (int c = 0; c < order + 1; c++) {
        A.set(0, c, 1);
      }
      double f = 0;
      for (int c = 0; c < order; c++) {
        f += deltas[c];
        double a = 1;
        for (int r = 0; r < order; r++) {
          a *= f / deltas[0];
          A.set(r + 1, c + 1, a);
        }
      }
      e.passEquationSys(&A, &x, &b);
      e.solve();
      for (int r = 0; r <= order; r++) {
        coefficients[r] = x.get(r);
      }
    }
    break;
  }
  case INTEGRATOR_ADAMSMOULTON: {
    // ADAMS-MOULTON order 1 to 6
    // right hand side vector
    for (int i = 0; i < order + 1; i++) {
      b.set(i, 1);
    }
    for (int i = 1; i < order + 1; i++) {
      A.set(i, 1, i); // second column
      A.set(1, i, 1); // second row
    }
    A.set(0, 0, 1);
    for (int c = 1; c <= order - 2; c++) {
      double entry = -c;
      for (int r = 2; r <= order; r++) {
        A.set(r, c + 2, r * entry);
        entry *= -c;
      }
    }
    e.passEquationSys(&A, &x, &b);
    e.solve();
    // vector x consists of a_{0}, b_{-1}, b_{0} ... b_{k-2} right here
    const double k = x.get(1);
    coefficients[0] = 1 / deltas[0] / k;
    coefficients[1] = -x.get(0) / deltas[0] / k;
    for (int i = 2; i <= order; i++) {
      coefficients[i] = -x.get(i) / k;
    }
    break;
  }
  }
}

/* Calculates the integration coefficient for numerical
   integration methods.  Supported methods are: Adams-Bashford (order
   1-6), forward Euler and explicit Gear (order 1-6). */
void transient::calcPredictorCoeff(const int method, const int order, const double *deltas,
                                   double *coefficients) {
  tmatrix<double> A(order + 1);
  tvector<double> x(order + 1);
  tvector<double> b(order + 1);
  eqnsys<double> e;
  e.setAlgo(ALGO_LU_DECOMPOSITION);

  switch (method) {
  case INTEGRATOR_EULER: {
    // FORWARD EULER
    coefficients[0] = 1;
    coefficients[1] = deltas[0];
    break;
  }
  case INTEGRATOR_GEAR: {
    // ADAMS-MOULTON order 1 to 6
    // right hand side vector
    b.set(0, 1);
    // first row
    for (int c = 0; c < order + 1; c++) {
      A.set(0, c, 1);
    }
    double f = 0;
    for (int c = 0; c < order + 1; c++) {
      f += deltas[c];
      double a = 1;
      for (int r = 0; r < order; r++) {
        a *= f / deltas[0];
        A.set(r + 1, c, a);
      }
    }
    e.passEquationSys(&A, &x, &b);
    e.solve();
    for (int r = 0; r <= order; r++) {
      coefficients[r] = x.get(r);
    }
    break;
  }
  case INTEGRATOR_ADAMSBASHFORD: {
    // ADAMS-BASHFORD order 1 to 6
    // right hand side vector
    for (int i = 0; i < order + 1; i++) {
      b.set(i, 1);
    }
    for (int i = 1; i < order + 1; i++) {
      A.set(1, i, 1); // second row
    }
    A.set(0, 0, 1);
    for (int c = 1; c <= order - 1; c++) {
      double entry = -c;
      for (int r = 2; r <= order; r++) {
        A.set(r, c + 1, r * entry);
        entry *= -c;
      }
    }
    e.passEquationSys(&A, &x, &b);
    e.solve();
    // vector x consists of a_{0}, b_{0}, b_{1} ... b_{k-1} right here
    coefficients[0] = x.get(0);
    for (int i = 1; i <= order; i++) {
      coefficients[i] = x.get(i) * deltas[0];
    }
    if constexpr (!FIXED_COEFF) {
      if (order == 2) {
        const double f = -deltas[0] / (2 * deltas[1]);
        coefficients[0] = 1;
        coefficients[1] = (1 - f) * deltas[0];
        coefficients[2] = f * deltas[0];
      }
    }
    break;
  }
  }
}

// Loads the equivalent conductance.
void transient::getConductance(integrator *c, double cap, double &geq) {
  const double *coeff = c->getCoefficients();
  geq = cap * coeff[0];
}

// Implicit Euler integrator.
void transient::integrateEuler(integrator *c, int qstate, double cap, double &geq, double &ceq) {
  const double *coeff = c->getCoefficients();
  int cstate = qstate + 1;
  geq = cap * coeff[0];
  ceq = c->getState(qstate, 1) * coeff[1];
  double cur = c->getState(qstate) * coeff[0] + ceq;
  c->setState(cstate, cur);
}

// Trapezoidal integrator.
void transient::integrateBilinear(integrator *c, int qstate, double cap, double &geq, double &ceq) {
  const double *coeff = c->getCoefficients();
  int cstate = qstate + 1;
  geq = cap * coeff[0];
  ceq = c->getState(qstate, 1) * coeff[1] - c->getState(cstate, 1);
  double cur = c->getState(qstate) * coeff[0] + ceq;
  c->setState(cstate, cur);
}

// Integrator using the Gear coefficients.
void transient::integrateGear(integrator *c, int qstate, double cap, double &geq, double &ceq) {
  const double *coeff = c->getCoefficients();
  int cstate = qstate + 1;
  geq = cap * coeff[0];
  ceq = 0;
  for (int i = 1; i <= c->getOrder(); i++) {
    ceq += c->getState(qstate, i) * coeff[i];
  }
  double cur = c->getState(qstate) * coeff[0] + ceq;
  c->setState(cstate, cur);
}

// Integrator using the Adams-Moulton coefficients.
void transient::integrateMoulton(integrator *c, int qstate, double cap, double &geq, double &ceq) {
  const double *coeff = c->getCoefficients();
  int cstate = qstate + 1;
  geq = cap * coeff[0];
  ceq = c->getState(qstate, 1) * coeff[1];
  for (int i = 2; i <= c->getOrder(); i++) {
    ceq += c->getState(cstate, i - 1) * coeff[i];
  }
  double cur = c->getState(qstate) * coeff[0] + ceq;
  c->setState(cstate, cur);
}

/* Applies the appropriate integration function to the given circuit object. */
void transient::setIntegrationMethod(integrator *c, int method) {
  switch (method) {
  case INTEGRATOR_GEAR:
    c->setIntegrateFunc(integrateGear);
    break;
  case INTEGRATOR_TRAPEZOIDAL:
    c->setIntegrateFunc(integrateBilinear);
    break;
  case INTEGRATOR_EULER:
    c->setIntegrateFunc(integrateEuler);
    break;
  case INTEGRATOR_ADAMSMOULTON:
    c->setIntegrateFunc(integrateMoulton);
    break;
  default:
    c->setIntegrateFunc(nullptr);
    break;
  }
  c->setConductorFunc(getConductance);
}

/* Returns an appropriate integrator type identifier and the maximum
   order depending on the given string argument. */
int transient::correctorType(const char *const method, int &maxOrder) {
  if (!strcmp(method, "Gear")) {
    if (maxOrder > 6)
      maxOrder = 6;
    if (maxOrder < 1)
      maxOrder = 1;
    return INTEGRATOR_GEAR;
  }
  if (!strcmp(method, "Trapezoidal")) {
    maxOrder = 2;
    return INTEGRATOR_TRAPEZOIDAL;
  }
  if (!strcmp(method, "Euler")) {
    maxOrder = 1;
    return INTEGRATOR_EULER;
  }
  if (!strcmp(method, "AdamsMoulton")) {
    if (maxOrder > 6)
      maxOrder = 6;
    if (maxOrder < 1)
      maxOrder = 1;
    return INTEGRATOR_ADAMSMOULTON;
  }
  if (!strcmp(method, "AdamsBashford")) {
    if (maxOrder > 6)
      maxOrder = 6;
    if (maxOrder < 1)
      maxOrder = 1;
    return INTEGRATOR_ADAMSBASHFORD;
  }
  return INTEGRATOR_UNKNOWN;
}

/* Returns the appropriate predictor integration method
   for the given corrector method and adjusts the order of the
   predictor as well based on the given corrector method. */
int transient::predictorType(int corrMethod, int corrOrder, int &predOrder) {
  int predMethod = INTEGRATOR_UNKNOWN;
  switch (corrMethod) {
  case INTEGRATOR_GEAR:
    predMethod = INTEGRATOR_GEAR;
    break;
  case INTEGRATOR_ADAMSMOULTON:
    predMethod = INTEGRATOR_ADAMSBASHFORD;
    break;
  case INTEGRATOR_TRAPEZOIDAL:
    predMethod = INTEGRATOR_ADAMSBASHFORD;
    break;
  case INTEGRATOR_EULER:
    predMethod = INTEGRATOR_EULER;
    break;
  }
  predOrder = corrOrder;
  return predMethod;
}

// Defines the integration algorithm for each possible order.
struct integration_types_t {
  int method;
  int integratorType[6];
  double corrErrorConstant[6];
  double predErrorConstant[6];
};

static integration_types_t integration_types[] = {
    {
        INTEGRATOR_EULER,
        {
            INTEGRATOR_EULER,
        },
        {
            -1.0 / 2,
        },
        {
            +1.0 / 2,
        },
    },
    {
        INTEGRATOR_TRAPEZOIDAL,
        {
            INTEGRATOR_EULER,
            INTEGRATOR_TRAPEZOIDAL,
        },
        {
            -1.0 / 2,
            -1.0 / 12,
        },
        {
            +1.0 / 2,
            +5.0 / 12,
        },
    },
    {
        INTEGRATOR_GEAR,
        {
            INTEGRATOR_GEAR,
            INTEGRATOR_GEAR,
            INTEGRATOR_GEAR,
            INTEGRATOR_GEAR,
            INTEGRATOR_GEAR,
            INTEGRATOR_GEAR,
        },
        {
            -1.0 / 2,
            -2.0 / 9,
            -3.0 / 22,
            -12.0 / 125,
            -10.0 / 137,
            -20.0 / 343,
        },
        {
            +1.0,
            +1.0,
            +1.0,
            +1.0,
            +1.0,
            +1.0,
        },
    },
    {
        INTEGRATOR_ADAMSMOULTON,
        {
            INTEGRATOR_ADAMSMOULTON,
            INTEGRATOR_ADAMSMOULTON,
            INTEGRATOR_ADAMSMOULTON,
            INTEGRATOR_ADAMSMOULTON,
            INTEGRATOR_ADAMSMOULTON,
            INTEGRATOR_ADAMSMOULTON,
        },
        {
            -1.0 / 2,
            -1.0 / 12,
            -1.0 / 24,
            -19.0 / 720,
            -3.0 / 160,
            -863.0 / 60480,
        },
        {
            +1.0 / 2,
            +1.0 / 12,
            +1.0 / 24,
            +19.0 / 720,
            +3.0 / 160,
            +863.0 / 60480,
        },
    },
    {
        INTEGRATOR_ADAMSBASHFORD,
        {
            INTEGRATOR_ADAMSBASHFORD,
            INTEGRATOR_ADAMSBASHFORD,
            INTEGRATOR_ADAMSBASHFORD,
            INTEGRATOR_ADAMSBASHFORD,
            INTEGRATOR_ADAMSBASHFORD,
            INTEGRATOR_ADAMSBASHFORD,
        },
        {
            -1.0 / 2,
            -5.0 / 12,
            -3.0 / 8,
            -251.0 / 720,
            -95.0 / 288,
            -19087.0 / 60480,
        },
        {
            +1.0 / 2,
            +5.0 / 12,
            +3.0 / 8,
            +251.0 / 720,
            +95.0 / 288,
            +19087.0 / 60480,
        },
    },
};

/* Returns the appropriate integration type for the given corrector integration type and order. */
int transient::correctorType(const int method, const int order) {
  return integration_types[method].integratorType[order - 1];
}

/* Returns the error constant for the given corrector. */
double transient::getCorrectorError(const int method, const int order) {
  return integration_types[method].corrErrorConstant[order - 1];
}

/* Returns the error constant for the given predictor. */
double transient::getPredictorError(const int method, const int order) {
  return integration_types[method].predErrorConstant[order - 1];
}

} // namespace qucs
