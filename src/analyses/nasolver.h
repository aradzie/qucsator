/*
 * nasolver.h - nodal analysis solver class definitions
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2008 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __NASOLVER_H__
#define __NASOLVER_H__

#include <string>

#include "analysis.h"
#include "eqnsys.h"
#include "tmatrix.h"
#include "tvector.h"

#define CONV_None 0
#define CONV_Attenuation 1
#define CONV_LineSearch 2
#define CONV_SteepestDescent 3
#define CONV_GMinStepping 4
#define CONV_SourceStepping 5

namespace qucs {

class analysis;
class circuit;
class nodelist;
class vector;

template <class nr_type_t> class nasolver : public analysis {
public:
  nasolver();
  explicit nasolver(const std::string &);
  nasolver(const nasolver &) = delete;
  ~nasolver() override;

  void setDescription(const std::string &n) { desc = n; }
  std::string getDescription() const { return desc; }

  void solve_pre();
  void solve_post();
  int solve_once();
  int solve_linear();
  int solve_nonlinear();
  int solve_nonlinear_continuation_gMinStepping();
  int solve_nonlinear_continuation_SourceStepping();

  void saveResults(const std::string &, const std::string &, int, qucs::vector *f = nullptr);

  typedef void (*calculate_func_t)(nasolver<nr_type_t> *);
  void setCalculation(calculate_func_t f) { calculate_func = f; }
  void calculate() {
    if (calculate_func) {
      (*calculate_func)(this);
    }
  }
  const char *getHelperDescription();

  // Returns the number of node voltages in the circuit.
  int getN();
  // Returns the number of branch currents in the circuit.
  int getM();

protected:
  /* Apply starting value to previous iteration value. */
  void restartDC();
  void savePreviousIteration();
  void restorePreviousIteration();
  int countNodes();
  int getNodeIndex(const std::string &);
  int findAssignedNode(circuit *, int);
  int countVoltageSources();
  circuit *findVoltageSource(int);
  void solveLinearEquations();
  bool checkConvergence();
  void saveSolution();

  void createMatrix();
  void createNoiseMatrix();

  void applyNodeset(bool reset = true);

private:
  void assignVoltageSources();

  void createGMatrix();
  void createBMatrix();
  void createCMatrix();
  void createDMatrix();
  void createIVector();
  void createEVector();
  void createZVector();

  void applyAttenuation();
  void lineSearch();
  void steepestDescent();

  std::string createV(int, const std::string &, int);
  std::string createI(int, const std::string &, int);
  std::string createOP(const std::string &, const std::string &);

  void saveNodeVoltages();
  void saveBranchCurrents();

  nr_type_t MatValX(nr_complex_t, nr_complex_t *);
  nr_type_t MatValX(nr_complex_t, double *);

protected:
  /* The right hand side input vector of the SLE. */
  tvector<nr_type_t> *z, *zprev;
  /* The solution vector of the SLE, the result of A*x=z. */
  tvector<nr_type_t> *x, *xprev;
  /* The left hand side input matrix of the SLE. */
  tmatrix<nr_type_t> *A;
  tmatrix<nr_type_t> *C;

  int iterations;
  int convHelper;
  int eqnAlgo;
  int updateMatrix;
  double gMin, srcFactor;
  std::string desc;
  nodelist *nlist; // ARA: This list exists for the duration of a single analysis.

private:
  eqnsys<nr_type_t> *eqns;
  double reltol;
  double abstol;
  double vntol;

private:
  calculate_func_t calculate_func;
};

} // namespace qucs

#include "nasolver.cpp"

#endif /* __NASOLVER_H__ */
