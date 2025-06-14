/*
 * trsolver.cpp - transient solver class implementation
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2009 Stefan Jahn <stefan@lkcc.org>
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

#include "trsolver.h"
#include "analysis.h"
#include "circuit.h"
#include "complex.h"
#include "constants.h"
#include "dataset.h"
#include "exception.h"
#include "exceptionstack.h"
#include "history.h"
#include "logging.h"
#include "nasolver.h"
#include "net.h"
#include "netdefs.h"
#include "object.h"
#include "sweep.h"
#include "transient.h"
#include "vector.h"

#define STEPDEBUG 0   // set to zero for release

namespace qucs {

using namespace transient;

trsolver::trsolver() {
  swp = nullptr;
  type = ANALYSIS_TRANSIENT;
  setDescription("transient");
  for (int i = 0; i < 8; i++) {
    solution[i] = nullptr;
  }
  hist = nullptr;
}

trsolver::trsolver(const std::string &n) : nasolver(n) {
  swp = nullptr;
  type = ANALYSIS_TRANSIENT;
  setDescription("transient");
  for (int i = 0; i < 8; i++) {
    solution[i] = nullptr;
  }
  hist = nullptr;
}

trsolver::~trsolver() {
  delete swp;
  for (int i = 0; i < 8; i++) {
    if (solution[i] != nullptr) {
      delete solution[i];
    }
  }
  delete hist;
}

int trsolver::solve() {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::solve()\n", getName());

  int error = 0, convError = 0;
  const char *const solver = getPropertyString("Solver");
  const bool initialDC = !strcmp(getPropertyString("initialDC"), "yes") ? true : false;

  runs++;
  double saveCurrent = current = 0;
  stepDelta = -1;
  converged = 0;
  statRejected = statSteps = statIterations = statConvergence = 0;

  if (!strcmp(solver, "CroutLU"))
    eqnAlgo = ALGO_LU_DECOMPOSITION;
  else if (!strcmp(solver, "DoolittleLU"))
    eqnAlgo = ALGO_LU_DECOMPOSITION_DOOLITTLE;
  else if (!strcmp(solver, "HouseholderQR"))
    eqnAlgo = ALGO_QR_DECOMPOSITION;
  else if (!strcmp(solver, "HouseholderLQ"))
    eqnAlgo = ALGO_QR_DECOMPOSITION_LS;
  else if (!strcmp(solver, "GolubSVD"))
    eqnAlgo = ALGO_SV_DECOMPOSITION;

  // Perform initial DC analysis.
  if (initialDC) {
    error = dcAnalysis();
    if (error) {
      return -1;
    }
  }

  // Initialize transient analysis.
  setDescription("transient");
  initTR();
  setCalculation((calculate_func_t)&calcTR);
  solve_pre();

  swp = createSweep("time");
  swp->reset();

  // Recall the DC solution.
  recallDcSolution();

  // Apply the nodesets and adjust previous solutions.
  applyNodeset(false);
  // Stores the given vector into all the solution vectors.
  for (int i = 0; i < 8; i++) {
    *solution[i] = *x;
  }

  // Tell integrators to be initialized.
  setMode(MODE_INIT);

  int running = 0;
  rejected = 0;
  delta /= 10;
  for (int i = 0; i < 8; i++) {
    deltas[i] = delta;
  }
  adjustOrder(1);

  // Start to sweep through time.
  for (int i = 0; i < swp->getSize(); i++) {
    const double time = swp->next();

#if DEBUG
    logprint(LOG_STATUS, "NOTIFY: %s: solving netlist for t = %e\n", getName(), time);
#endif

    do // while (saveCurrent < time), i.e. until a requested breakpoint is hit
    {
      logprint(LOG_STATUS, "DEBUG: %s: t = %.3e\n", getName(), current);

#if STEPDEBUG
      if (delta == deltaMin) {
        // the integrator step size has become smaller than the
        // specified allowed minimum, Qucs is unable to solve the circuit
        // while meeting the tolerance conditions
        logprint(LOG_ERROR, "WARNING: %s: minimum delta h = %.3e at t = %.3e\n", getName(), delta,
                 current);
      }
#endif
      // updates the integrator coefficients, and updates the array of prev
      // 8 deltas with the new delta for this step
      deltas[0] = delta;
      calcCorrectorCoeff(corrMethod, corrOrder, deltas, corrCoeff);
      calcPredictorCoeff(predMethod, predOrder, deltas, predCoeff);

      // Run predictor to get a start value for the solution vector for
      // the successive iterative corrector process
      predictor();

      // restart Newton iteration
      if (rejected) {
        restartDC(); // restart non-linear devices
        rejected = 0;
      }

      // Run corrector process with appropriate exception handling.
      // The corrector iterates through the solutions of the integration
      // process until a certain error tolerance has been reached.
      error = corrector();

      if (estack.top()) {
        switch (estack.top()->getCode()) {
        case EXCEPTION_NO_CONVERGENCE:
          estack.pop();

          // step back from the current time value to the previous time
          if (current > 0) {
            current -= delta;
          }
          // Reduce step-size (by half) if failed to converge.
          delta /= 2;
          if (delta <= deltaMin) {
            // but do not reduce the step size below a specified minimum
            delta = deltaMin;
            // instead reduce the order of the integration
            adjustOrder(1);
          }
          // step forward to the new current time value
          if (current > 0)
            current += delta;

          // Update statistics.
          statRejected++;
          statConvergence++;
          rejected++; // mark the previous step size choice as rejected
          converged = 0;
          error = 0;

          // Start using damped Newton-Raphson.
          convHelper = CONV_SteepestDescent;
          convError = 2;
          break;
        default:
          // Otherwise return.
          estack.print();
          return -1;
        }
      }

      // return if any errors occurred other than convergence failure
      if (error) {
        return -1;
      }

      // if the step was rejected, the solution loop is restarted here
      if (rejected) {
        continue;
      }

      // check whether Jacobian matrix is still non-singular
      if (!A->isFinite()) {
        logprint(LOG_ERROR,
                 "ERROR: %s: Jacobian singular at t = %.3e, "
                 "aborting %s analysis\n",
                 getName(), current, getDescription().c_str());
        return -1;
      }

      // Update statistics and no more damped Newton-Raphson.
      statIterations += iterations;
      if (--convError < 0) {
        convHelper = 0;
      }

      // Now advance in time or not...
      if (running > 1) {
        adjustDelta(time); // Will call `nextStates()`.
        adjustOrder();
      } else {
        fillStates();
        nextStates();
        rejected = 0;
      }

      saveCurrent = current;
      current += delta;
      running++;
      converged++;

      // Tell integrators to be running.
      setMode(MODE_NONE);

      // Initialize or update history.
      if (running > 1) {
        updateHistory(saveCurrent);
      } else {
        initHistory(saveCurrent);
      }
    } while (saveCurrent < time);

    // Save results.
#if STEPDEBUG
    logprint(LOG_STATUS, "DEBUG: save point at t = %.3e, h = %.3e\n", saveCurrent, delta);
#endif

    saveAllResults(time);
  } // for (int i = 0; i < swp->getSize (); i++)

  solve_post();

  logprint(LOG_STATUS, "NOTIFY: %s: average time-step %g, %d rejections\n", getName(),
           (saveCurrent / statSteps), statRejected);
  logprint(LOG_STATUS,
           "NOTIFY: %s: average NR-iterations %g, "
           "%d non-convergences\n",
           getName(), statIterations / statSteps, statConvergence);

  return 0;
}

/* Predicts a start value for the solution vector for the successive iterative corrector process. */
void trsolver::predictor() {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::predictor()\n", getName());
  switch (predMethod) {
  case INTEGRATOR_EULER: // explicit forward EULER
    predictEuler();
    break;
  case INTEGRATOR_GEAR: // explicit GEAR
    predictGear();
    break;
  case INTEGRATOR_ADAMSBASHFORD: // ADAMS-BASHFORD
    predictBashford();
    break;
  default:
    *x = *solution[1]; // This too is a simple predictor...
    break;
  }
  saveSolution();
  *solution[0] = *x;
}

/* Predicts the successive solution vector using the explicit forward Euler integration formula.
 * Actually this is Adams-Bashford order 1. */
void trsolver::predictEuler() {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::predictEuler()\n", getName());
  const int N = countNodes();
  const int M = countVoltageSources();
  // go through each solution
  for (int r = 0; r < N + M; r++) {
    double xn = predCoeff[0] * solution[1]->get(r);
    const double hn = deltas[1];
    const double dd = (solution[1]->get(r) - solution[2]->get(r)) / hn;
    xn += predCoeff[1] * dd;
    x->set(r, xn); // save prediction
  }
}

/* Predicts the successive solution vector using the explicit Gear integration formula. */
void trsolver::predictGear() {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::predictGear()\n", getName());
  const int N = countNodes();
  const int M = countVoltageSources();
  // go through each solution
  for (int r = 0; r < N + M; r++) {
    double xn = 0;
    for (int o = 0; o <= predOrder; o++) {
      // a0, a1, ... coefficients
      xn += predCoeff[o] * solution[o + 1]->get(r);
    }
    x->set(r, xn); // save prediction
  }
}

/* Predicts the successive solution vector using the explicit Adams-Bashford integration formula. */
void trsolver::predictBashford() {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::predictBashford()\n", getName());
  const int N = countNodes();
  const int M = countVoltageSources();
  // go through each solution
  for (int r = 0; r < N + M; r++) {
    double xn = predCoeff[0] * solution[1]->get(r); // a0 coefficient
    for (int o = 1; o <= predOrder; o++) {
      const double hn = deltas[o]; // previous time-step
      // divided differences
      const double dd = (solution[o]->get(r) - solution[o + 1]->get(r)) / hn;
      xn += predCoeff[o] * dd; // b0, b1, ... coefficients
    }
    x->set(r, xn); // save prediction
  }
}

/* Iterates through the solutions of the integration process
 * until a certain error tolerance has been reached. */
int trsolver::corrector() {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::corrector()\n", getName());

  return solve_nonlinear();
}

// The function advances one more time-step.
void trsolver::nextStates() {
  circuit *root = subnet->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    c->nextState();
  }
  // Save current solution.
  *solution[0] = *x;
  // Shift solutions, free space for the next solution vector X.
  for (int i = 7; i > 0; i--) {
    *solution[i] = *solution[i - 1];
    deltas[i] = deltas[i - 1];
  }
  statSteps++;
}

/* Stores the current state of each circuit into all other states as well.
 * It is useful for higher order integration methods in order to initialize the states
 * after the initial transient solution. */
void trsolver::fillStates() {
  circuit *root = subnet->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    for (int s = 0; s < c->getStates(); s++) {
      c->fillState(s, c->getState(s));
    }
  }
}

// Modifies the circuit lists integrator mode.
void trsolver::setMode(int mode) {
  circuit *root = subnet->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    c->setMode(mode);
  }
}

/* Adapts the current time-step the transient analysis.
 * For the computation of the new time-step the truncation error depending
 * on the integration method is used. */
double trsolver::checkDelta() {
  double LTEreltol = getPropertyDouble("LTEreltol");
  double LTEabstol = getPropertyDouble("LTEabstol");
  double LTEfactor = getPropertyDouble("LTEfactor");
  double dif, rel, tol, lte, q, n = std::numeric_limits<double>::max();
  const int N = countNodes();
  const int M = countVoltageSources();

  // cec = corrector error constant
  double cec = getCorrectorError(corrMethod, corrOrder);
  // pec = predictor error constant
  double pec = getPredictorError(predMethod, predOrder);

  // go through each solution
  for (int r = 0; r < N + M; r++) {

    // skip real voltage sources
    if (r >= N) {
      if (findVoltageSource(r - N)->isVSource())
        continue;
    }

    dif = x->get(r) - solution[0]->get(r);
    if (std::isfinite(dif) && dif != 0) {
      // use Milne' estimate for the local truncation error
      rel = MAX(fabs(x->get(r)), fabs(solution[0]->get(r)));
      tol = LTEreltol * rel + LTEabstol;
      lte = LTEfactor * (cec / (pec - cec)) * dif;
      q = delta * exp(log(fabs(tol / lte)) / (corrOrder + 1));
      n = std::min(n, q);
    }
  }
#if STEPDEBUG
  logprint(LOG_STATUS,
           "DEBUG: delta according to local truncation "
           "error h = %.3e\n",
           n);
#endif
  delta = std::min((n > 1.9 * delta) ? 2 * delta : delta, n);
  return delta;
}

/* Adapts the current time-step according to the global truncation error. */
void trsolver::adjustDelta(double t) {
  const bool relaxTSR = !strcmp(getPropertyString("relaxTSR"), "yes") ? true : false;

  deltaOld = delta;
  delta = checkDelta();
  if (delta > deltaMax) {
    delta = deltaMax;
  }
  if (delta < deltaMin) {
    delta = deltaMin;
  }

  // delta correction in order to hit exact breakpoint
  int good = 0;
  // relaxed step raster?
  if (!relaxTSR) {
    // Is this a good guess?
    if (!statConvergence || converged > 64) {
      // check next breakpoint
      if (stepDelta > 0.0) {
        // restore last valid delta
        delta = stepDelta;
        stepDelta = -1.0;
      } else {
        // check whether this step will bring too close to a breakpoint
        //   is ok if the step will go past the breakpoint, this is handled by
        //   the next branch
        if (t - (current + delta) < deltaMin && current + delta < t) {
          // if we take this delta we will end up too close to the breakpoint
          // and next step will be very tiny, possibly causing numerical issues
          // so reduce it so that next step will likely not end up too close
          // to the breakpoint
          delta /= 2.0;
        } else {
          if (delta > (t - current) && t > current) {
            // save last valid delta and set exact step
            stepDelta = deltaOld;
            delta = t - current;
            good = 1;
          } else {
            stepDelta = -1.0;
          }
        }
      }
      if (delta > deltaMax)
        delta = deltaMax;
      if (delta < deltaMin)
        delta = deltaMin;
    }
  }

  // usual delta correction
  if (delta > 0.9 * deltaOld || good) {
    // accept current delta
    nextStates();
    rejected = 0;
#if STEPDEBUG
    logprint(LOG_STATUS, "DEBUG: delta accepted at t = %.3e, h = %.3e\n", current, delta);
#endif
  } else if (deltaOld > delta) {
    // reject current delta
    rejected++;
    statRejected++;
#if STEPDEBUG
    logprint(LOG_STATUS, "DEBUG: delta rejected at t = %.3e, h = %.3e\n", current, delta);
#endif
    if (current > 0)
      current -= deltaOld;
  } else {
    nextStates();
    rejected = 0;
  }
}

/* Increase or reduces the current order of the integration method. */
void trsolver::adjustOrder(int reduce) {
  if ((corrOrder < corrMaxOrder && !rejected) || reduce) {
    if (reduce) {
      corrOrder = 1;
    } else if (!rejected) {
      corrOrder++;
    }

    // adjust type and order of corrector and predictor
    corrMethod = correctorType(corrMethod0, corrOrder);
    predMethod = predictorType(corrMethod, corrOrder, predOrder);

    // apply new corrector method and order to each circuit
    circuit *root = subnet->getRoot();
    for (circuit *c = root; c != nullptr; c = c->getNext()) {
      c->setOrder(corrOrder);
      setIntegrationMethod(c, corrMethod);
    }
  }
}

/* Goes through the list of circuit objects and runs its initTR() function. */
void trsolver::initTR() {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::initTR()\n", getName());

  const char *const IMethod = getPropertyString("IntegrationMethod");
  double start = getPropertyDouble("Start");
  double stop = getPropertyDouble("Stop");
  double points = getPropertyDouble("Points");

  // fetch corrector integration method and determine predicor method
  corrMaxOrder = getPropertyInteger("Order");
  corrMethod = corrMethod0 = correctorType(IMethod, corrMaxOrder);
  predMethod = predictorType(corrMethod, corrMaxOrder, predMaxOrder);
  corrOrder = corrMaxOrder;
  predOrder = predMaxOrder;

  // initialize step values
  delta = getPropertyDouble("InitialStep");
  deltaMin = getPropertyDouble("MinStep");
  deltaMax = getPropertyDouble("MaxStep");
  if (deltaMax == 0.0) {
    deltaMax = std::min((stop - start) / (points - 1), stop / 200);
  }
  if (deltaMin == 0.0) {
    deltaMin = NR_TINY * 10 * deltaMax;
  }
  if (delta == 0.0) {
    delta = std::min(stop / 200, deltaMax) / 10;
  }
  if (delta < deltaMin) {
    delta = deltaMin;
  }
  if (delta > deltaMax) {
    delta = deltaMax;
  }

  // initialize history of solution vectors (solutions)
  for (int i = 0; i < 8; i++) {
    // solution contains the last sets of node voltages and branch
    // currents at each of the last 8 'deltas'.
    solution[i] = new tvector<double>;
    deltas[i] = delta;
  }

  // set the initial corrector and predictor coefficients
  calcCorrectorCoeff(corrMethod, corrOrder, deltas, corrCoeff);
  calcPredictorCoeff(predMethod, predOrder, deltas, predCoeff);

  // tell circuits about the transient analysis
  circuit *c, *root = subnet->getRoot();
  for (c = root; c != nullptr; c = c->getNext()) {
    c->initTR();
    // Share the state with the circuit.
    c->setDelta(deltas); // ARA: Only used by the BJT model.
    c->setCoefficients(corrCoeff);
    c->setOrder(corrOrder);
    setIntegrationMethod(c, corrMethod);
  }
  // also initialize created circuits
  for (c = root; c != nullptr; c = c->getPrev()) {
    c->initTR();
    // Share the state with the circuit.
    c->setDelta(deltas); // ARA: Only used by the BJT model.
    c->setCoefficients(corrCoeff);
    c->setOrder(corrOrder);
    setIntegrationMethod(c, corrMethod);
  }
}

/* Goes through the list of circuit objects and runs its calcTR() function. */
void trsolver::calcTR(trsolver *self) {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::calcTR()\n", self->getName());

  circuit *root = self->getNet()->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    c->calcTR(self->current);
  }
}

void trsolver::initHistory(double t) {
  // initialize time vector
  hist = new history();
  hist->push_back(t);
  hist->self();
  // initialize circuit histories
  double age = 0.0;
  circuit *root = subnet->getRoot();
  for (circuit *c = root; c != nullptr; c = c->getNext()) {
    if (c->hasHistory()) {
      c->applyHistory(hist);
      saveHistory(c);
      if (c->getHistoryAge() > age) {
        age = c->getHistoryAge();
      }
    }
  }
  // set maximum required age for all circuits
  hist->setAge(age);
}

/* Updates the histories for the circuits which requested them. */
void trsolver::updateHistory(double t) {
  if (t > hist->last()) {
    // update time vector
    hist->push_back(t);
    // update circuit histories
    circuit *root = subnet->getRoot();
    for (circuit *c = root; c != nullptr; c = c->getNext()) {
      if (c->hasHistory())
        saveHistory(c);
    }
    hist->drop();
  }
}

// Stores node voltages and branch currents in the given circuits history.
void trsolver::saveHistory(circuit *c) {
  const int N = countNodes();
  const int size = c->getSize();
  for (int i = 0; i < size; i++) {
    // save node voltages
    if (int r = findAssignedNode(c, i); r < 0) {
      // the node was not found, append a zero to the history
      // matching this index
      c->appendHistory(i, 0.0);
    } else {
      // the node was found, append the voltage value to
      // that node's history
      c->appendHistory(i, x->get(r));
    }
  }

  for (int i = 0; i < c->getVoltageSources(); i++) {
    // save branch currents
    const int r = c->getVoltageSource() + i;
    c->appendHistory(i + size, x->get(r + N));
  }
}

/* Saves the results of a single solve() functionality
   (for the given timestamp) into the output dataset. */
void trsolver::saveAllResults(double time) {
  logprint(LOG_STATUS, "NOTIFY: %s: trsolver::saveAllResults(%e)\n", getName(), time);

  qucs::vector *t = data->findDependency("time");
  // add current frequency to the dependency of the output dataset
  if (t == nullptr) {
    data->addDependency(t = new qucs::vector("time"));
  }
  if (runs == 1) {
    t->add(time);
  }
  saveResults("Vt", "It", 0, t);
}

PROP_REQ[] = {
    {"Type", PROP_STR, {PROP_NO_VAL, "lin"}, PROP_RNG_STR2("lin", "log")},
    {"Start", PROP_REAL, {0, PROP_NO_STR}, PROP_POS_RANGE},
    {"Stop", PROP_REAL, {1e-3, PROP_NO_STR}, PROP_POS_RANGE},
    {"Points", PROP_INT, {10, PROP_NO_STR}, PROP_MIN_VAL(2)},
    PROP_NO_PROP,
};
PROP_OPT[] = {
    {"IntegrationMethod",
     PROP_STR,
     {PROP_NO_VAL, "Trapezoidal"},
     PROP_RNG_STR4("Euler", "Trapezoidal", "Gear", "AdamsMoulton")},
    {"Order", PROP_INT, {2, PROP_NO_STR}, PROP_RNGII(1, 6)},
    {"InitialStep", PROP_REAL, {1e-9, PROP_NO_STR}, PROP_POS_RANGE},
    {"MinStep", PROP_REAL, {1e-16, PROP_NO_STR}, PROP_POS_RANGE},
    {"MaxStep", PROP_REAL, {0, PROP_NO_STR}, PROP_POS_RANGE},
    {"MaxIter", PROP_INT, {150, PROP_NO_STR}, PROP_RNGII(2, 10000)},
    {"abstol", PROP_REAL, {1e-12, PROP_NO_STR}, PROP_RNG_X01I},
    {"vntol", PROP_REAL, {1e-6, PROP_NO_STR}, PROP_RNG_X01I},
    {"reltol", PROP_REAL, {1e-3, PROP_NO_STR}, PROP_RNG_X01I},
    {"LTEabstol", PROP_REAL, {1e-6, PROP_NO_STR}, PROP_RNG_X01I},
    {"LTEreltol", PROP_REAL, {1e-3, PROP_NO_STR}, PROP_RNG_X01I},
    {"LTEfactor", PROP_REAL, {1, PROP_NO_STR}, PROP_RNGII(1, 16)},
    {"Temp", PROP_REAL, {26.85, PROP_NO_STR}, PROP_MIN_VAL(K)},
    {"Solver", PROP_STR, {PROP_NO_VAL, "CroutLU"}, PROP_RNG_SOL},
    {"relaxTSR", PROP_STR, {PROP_NO_VAL, "no"}, PROP_RNG_YESNO},
    {"initialDC", PROP_STR, {PROP_NO_VAL, "yes"}, PROP_RNG_YESNO},
    PROP_NO_PROP,
};
struct define_t trsolver::anadef = {"TR", 0, PROP_ACTION, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF};

} // namespace qucs
