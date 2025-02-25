/*
 * circuit.cpp - circuit class implementation
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
 */

#include <cassert>
#include <cstring>

#include "circuit.h"
#include "complex.h"
#include "component_id.h"
#include "history.h"
#include "logging.h"
#include "matrix.h"
#include "microstrip/substrate.h"
#include "node.h"
#include "object.h"
#include "tvector.h"
#include "valuelist.h"

namespace qucs {

// normalising impedance
const double circuit::z0 = 50.0;

circuit::circuit() : object(), integrator() {
  next = prev = nullptr;
  size = 0;
  MatrixN = MatrixS = MatrixY = nullptr;
  MatrixB = MatrixC = MatrixD = nullptr;
  VectorQ = VectorE = VectorI = VectorV = VectorJ = nullptr;
  MatrixQV = nullptr;
  VectorCV = VectorGV = nullptr;
  nodes = nullptr;
  pacport = 0;
  pol = 1;
  flag = CIRCUIT_ORIGINAL | CIRCUIT_LINEAR;
  subst = nullptr;
  vsource = -1;
  vsources = 0;
  nsources = 0;
  inserted = -1;
  subcircuit = std::string();
  subnet = nullptr;
  deltas = nullptr;
  histories = nullptr;
  nHistories = 0;
  type = CIR_UNKNOWN;
}

circuit::circuit(int s) : object(), integrator() {
  next = prev = nullptr;
  assert(s >= 0);
  size = s;
  if (size > 0)
    nodes = new node[s];
  MatrixN = MatrixS = MatrixY = nullptr;
  MatrixB = MatrixC = MatrixD = nullptr;
  VectorQ = VectorE = VectorI = VectorV = VectorJ = nullptr;
  MatrixQV = nullptr;
  VectorCV = VectorGV = nullptr;
  pacport = 0;
  pol = 1;
  flag = CIRCUIT_ORIGINAL | CIRCUIT_LINEAR;
  subst = nullptr;
  vsource = -1;
  vsources = 0;
  nsources = 0;
  inserted = -1;
  subcircuit = std::string();
  subnet = nullptr;
  deltas = nullptr;
  histories = nullptr;
  nHistories = 0;
  type = CIR_UNKNOWN;
}

circuit::~circuit() {
  if (size > 0) {
    delete[] MatrixS;
    delete[] MatrixN;
    freeMatrixMNA();
    freeMatrixHB();
    delete[] nodes;
  }
  deleteHistory();
}

/* With this function the number of ports of the circuit object can be
   changed.  Previously stored node and matrix information gets
   completely lost except the current size equals the given size. */
void circuit::setSize(int s) {
  // nothing to do here
  if (size == s)
    return;
  assert(s >= 0);

  if (size > 0) {
    // destroy any matrix and node information
    delete[] MatrixS;
    delete[] MatrixN;
    MatrixS = MatrixN = nullptr;
    freeMatrixMNA();
    delete[] nodes;
    nodes = nullptr;
  }

  if ((size = s) > 0) {
    // re-create matrix and node information space
    nodes = new node[size];
    allocMatrixS();
    allocMatrixN(nsources);
    allocMatrixMNA();
  }
}

/* Destroys the HB-matrix memory. */
void circuit::freeMatrixHB() {
  if (VectorQ) {
    delete[] VectorQ;
    VectorQ = nullptr;
  }
  if (MatrixQV) {
    delete[] MatrixQV;
    MatrixQV = nullptr;
  }
  if (VectorCV) {
    delete[] VectorCV;
    VectorCV = nullptr;
  }
  if (VectorGV) {
    delete[] VectorGV;
    VectorGV = nullptr;
  }
}

/* Allocates the HB-matrix memory. */
void circuit::allocMatrixHB() {
  if (VectorQ) {
    memset(VectorQ, 0, size * sizeof(nr_complex_t));
  } else {
    VectorQ = new nr_complex_t[size];
  }
  if (MatrixQV) {
    memset(MatrixQV, 0, size * size * sizeof(nr_complex_t));
  } else {
    MatrixQV = new nr_complex_t[size * size];
  }
  if (VectorCV) {
    memset(VectorCV, 0, size * sizeof(nr_complex_t));
  } else {
    VectorCV = new nr_complex_t[size];
  }
  if (VectorGV) {
    memset(VectorGV, 0, size * sizeof(nr_complex_t));
  } else {
    VectorGV = new nr_complex_t[size];
  }
}

/* Allocates the S-parameter matrix memory. */
void circuit::allocMatrixS() {
  if (MatrixS) {
    memset(MatrixS, 0, size * size * sizeof(nr_complex_t));
  } else {
    MatrixS = new nr_complex_t[size * size];
  }
}

/* Allocates the noise correlation matrix memory. */
void circuit::allocMatrixN(int sources) {
  nsources = sources;
  delete[] MatrixN;
  MatrixN = new nr_complex_t[(size + sources) * (size + sources)];
}

/* Allocates the matrix memory for the MNA matrices. */
void circuit::allocMatrixMNA() {
  freeMatrixMNA();
  if (size > 0) {
    MatrixY = new nr_complex_t[size * size]; // ARA: Input matrix, conductances.
    VectorI = new nr_complex_t[size]; // ARA: Input vector, currents.
    VectorV = new nr_complex_t[size]; // ARA: Output vector, voltages.
    if (vsources > 0) {
      MatrixB = new nr_complex_t[vsources * size]; // ARA: Input matrix.
      MatrixC = new nr_complex_t[vsources * size]; // ARA: Input matrix.
      MatrixD = new nr_complex_t[vsources * vsources]; // ARA: Input matrix.
      VectorE = new nr_complex_t[vsources]; // ARA: Input vector, voltages.
      VectorJ = new nr_complex_t[vsources]; // ARA: Output vector, currents.
    }
  }
}

/* Free()'s all memory used by the MNA matrices. */
void circuit::freeMatrixMNA() {
  if (MatrixY) {
    delete[] MatrixY;
    MatrixY = nullptr;
  }
  if (MatrixB) {
    delete[] MatrixB;
    MatrixB = nullptr;
  }
  if (MatrixC) {
    delete[] MatrixC;
    MatrixC = nullptr;
  }
  if (MatrixD) {
    delete[] MatrixD;
    MatrixD = nullptr;
  }
  if (VectorE) {
    delete[] VectorE;
    VectorE = nullptr;
  }
  if (VectorI) {
    delete[] VectorI;
    VectorI = nullptr;
  }
  if (VectorV) {
    delete[] VectorV;
    VectorV = nullptr;
  }
  if (VectorJ) {
    delete[] VectorJ;
    VectorJ = nullptr;
  }
}

/* Sets the name and port number of one of the circuit's
   nodes.  It also tells the appropriate node about the circuit it
   belongs to.  The optional 'intern' argument is used to mark a node
   to be for internal use only. */
void circuit::setNode(int i, const std::string &n, int intern) {
  nodes[i].setName(n);
  nodes[i].setCircuit(this);
  nodes[i].setPort(i);
  nodes[i].setInternal(intern);
}

node *circuit::getNode(int i) const { return &nodes[i]; }

#if DEBUG
// DEBUG function:  Prints the S parameters of this circuit object.
void circuit::print() {
  for (int i = 0; i < getSize(); i++) {
    for (int j = 0; j < getSize(); j++) {
      logprint(LOG_STATUS, "%s S%d%d(%+.3e,%+.3e) ", getName(), i, j, (double)real(getS(i, j)),
               (double)imag(getS(i, j)));
    }
    logprint(LOG_STATUS, "\n");
  }
}
#endif /* DEBUG */

/* Adds a operating point consisting of a key and a value to the circuit. */
void circuit::addOperatingPoint(const std::string &n, double val) {
  qucs::pair p(n, val);
  oper.insert({{n, p}});
}

/* Returns the requested operating point value which has been
   previously added as its double representation.  If there is no such
   operating point the function returns zero. */
double circuit::getOperatingPoint(const std::string &n) {
  const auto it = oper.find(n);
  if (it != oper.end())
    return (*it).second.getValue();
  return 0.0;
}

/* This function sets the operating point specified by the given name
   to the value passed to the function. */
void circuit::setOperatingPoint(const std::string &n, double val) {
  auto it = oper.find(n);
  if (it != oper.end())
    (*it).second.setValue(val);
  else
    addOperatingPoint(n, val);
}

/* The function checks whether the circuit has got a certain operating
   point value.  If so it returns non-zero, otherwise it returns
   zero. */
int circuit::hasOperatingPoint(const std::string &n) { return oper.find(n) != oper.end(); }

/* This function adds a characteristic point consisting of a key and a
   value to the circuit. */
void circuit::addCharacteristic(const std::string &n, double val) {
  qucs::pair p(n, val);
  charac.insert({{n, p}});
}

/* Returns the requested characteristic value which has been
   previously added as its double representation.  If there is no such
   characteristic value the function returns zero. */
double circuit::getCharacteristic(const std::string &n) {
  const auto it = charac.find(n);
  if (it != charac.end())
    return (*it).second.getValue();
  return 0.0;
}

/* This function sets the characteristic value specified by the given
   name to the value passed to the function. */
void circuit::setCharacteristic(const std::string &n, double val) {
  auto it = charac.find(n);
  if (it != charac.end())
    (*it).second.setValue(val);
  else
    addCharacteristic(n, val);
}

/* The function checks whether the circuit has got a certain
   characteristic value.  If so it returns non-zero, otherwise it
   returns zero. */
int circuit::hasCharacteristic(const std::string &n) { return charac.find(n) != charac.end(); }

/* The function returns an internal node or circuit name with the
   given prefix and based on the given circuits name.  The caller is
   responsible to free() the returned string. */
std::string circuit::createInternal(const std::string &prefix, const std::string &obj) {
  return "_" + prefix + "#" + obj;
}

/* Creates an internal node given the node number as well as the name
   suffix.  An appropriate node name is constructed from the circuits
   name and the suffix. */
void circuit::setInternalNode(int node, const std::string &suffix) {
  const std::string &n = createInternal(getName(), suffix);
  setNode(node, n, 1);
}

/* This function copies the matrix elements inside the given matrix to
   the internal S-parameter matrix of the circuit. */
void circuit::setMatrixS(matrix s) {
  int r = s.getRows();
  int c = s.getCols();
  // copy matrix elements
  if (r > 0 && c > 0 && r * c == size * size) {
    memcpy(MatrixS, s.getData(), sizeof(nr_complex_t) * r * c);
  }
}

/* The function return a matrix containing the S-parameters of the
   circuit. */
matrix circuit::getMatrixS() {
  matrix res(size);
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      res(i, j) = MatrixS[i * size + j];
  return res;
}

/* This function copies the matrix elements inside the given matrix to
   the internal noise correlation matrix of the circuit. */
void circuit::setMatrixN(matrix n) {
  int r = n.getRows();
  int c = n.getCols();
  // copy matrix elements
  if (r > 0 && c > 0 && r * c == size * size) {
    memcpy(MatrixN, n.getData(), sizeof(nr_complex_t) * r * c);
  }
}

/* The function return a matrix containing the noise correlation
   matrix of the circuit. */
matrix circuit::getMatrixN() {
  matrix res(size);
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      res(i, j) = MatrixN[i * size + j];
  return res;
}

/* This function copies the matrix elements inside the given matrix to
   the internal G-MNA matrix of the circuit. */
void circuit::setMatrixY(matrix y) {
  int r = y.getRows();
  int c = y.getCols();
  // copy matrix elements
  if (r > 0 && c > 0 && r * c == size * size) {
    memcpy(MatrixY, y.getData(), sizeof(nr_complex_t) * r * c);
  }
}

/* The function return a matrix containing the G-MNA matrix of the
   circuit. */
matrix circuit::getMatrixY() {
  matrix res(size);
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      res(i, j) = MatrixY[i * size + j];
  return res;
}

/* This function can be used by several components in order to place
   the n-th voltage source between node 'pos' and node 'neg' with the
   given value.  Remember to indicate this voltage source using the
   function setVoltageSources(). */
void circuit::voltageSource(int n, int pos, int neg, double value) {
  setC(n, pos, +1.0);
  setC(n, neg, -1.0);
  setB(pos, n, +1.0);
  setB(neg, n, -1.0);
  setD(n, n, 0.0);
  setE(n, value);
}

/* Runs the necessary calculation in order to perform a single integration step
 * of a voltage controlled capacitance placed in between the given nodes.
 * It is assumed that the appropriate charge only depends on the voltage between these nodes. */
void circuit::transientCapacitance(int qstate, int pos, int neg, double cap, double voltage,
                                   double charge) {
  double g, i;
  int cstate = qstate + 1;
  setState(qstate, charge);
  integrate(qstate, cap, g, i);
  addY(pos, pos, +g);
  addY(neg, neg, +g);
  addY(pos, neg, -g);
  addY(neg, pos, -g);
  i = pol * (getState(cstate) - g * voltage);
  addI(pos, -i);
  addI(neg, +i);
}

/* This is the one-node variant of the above function.  It performs
   the same steps for a single node related to ground. */
void circuit::transientCapacitance(int qstate, int node, double cap, double voltage,
                                   double charge) {
  double g, i;
  int cstate = qstate + 1;
  setState(qstate, charge);
  integrate(qstate, cap, g, i);
  addY(node, node, +g);
  i = pol * (getState(cstate) - g * voltage);
  addI(node, -i);
}

/* The function performs a single integration step of the given charge
   located between the given nodes.  It saves the current
   contributions of the charge itself and considers the polarity of
   the circuit. */
void circuit::transientCapacitanceQ(int qstate, int qpos, int qneg, double charge) {
  double unused, i;
  int cstate = qstate + 1;
  setState(qstate, charge);
  integrate(qstate, 0, unused, unused);
  i = pol * getState(cstate);
  addI(qpos, -i);
  addI(qneg, +i);
}

/* This is the one-node variant of the above function.  It performs
   the same steps for a single node related to ground. */
void circuit::transientCapacitanceQ(int qstate, int qpos, double charge) {
  double unused, i;
  int cstate = qstate + 1;
  setState(qstate, charge);
  integrate(qstate, 0, unused, unused);
  i = pol * getState(cstate);
  addI(qpos, -i);
}

/* This function stores the Jacobian entries due to the C = dQ/dV
   value.  The nodes where the charge is located as well as those of
   the voltage dependency, the appropriate capacitance value and the
   voltage across the the controlling branch must be given.  It also
   saves the current contributions which are necessary for the NR
   iteration and considers the polarity of the circuit. */
void circuit::transientCapacitanceC(int qpos, int qneg, int vpos, int vneg, double cap,
                                    double voltage) {
  double g, i;
  conductor(cap, g);
  addY(qpos, vpos, +g);
  addY(qneg, vneg, +g);
  addY(qpos, vneg, -g);
  addY(qneg, vpos, -g);
  i = pol * (g * voltage);
  addI(qpos, +i);
  addI(qneg, -i);
}

/* This is the one-node variant of the transientCapacitanceC()
   function.  It performs the same steps for a single charge node
   related to ground. */
void circuit::transientCapacitanceC2V(int qpos, int vpos, int vneg, double cap, double voltage) {
  double g, i;
  conductor(cap, g);
  addY(qpos, vpos, +g);
  addY(qpos, vneg, -g);
  i = pol * (g * voltage);
  addI(qpos, +i);
}

/* This is the one-node variant of the transientCapacitanceC()
   function.  It performs the same steps for a single voltage node
   related to ground. */
void circuit::transientCapacitanceC2Q(int qpos, int qneg, int vpos, double cap, double voltage) {
  double g, i;
  conductor(cap, g);
  addY(qpos, vpos, +g);
  addY(qneg, vpos, -g);
  i = pol * (g * voltage);
  addI(qpos, +i);
  addI(qneg, -i);
}

/* This is the one-node variant of the transientCapacitanceC()
   function.  It performs the same steps for a single voltage node and
   charge node related to ground. */
void circuit::transientCapacitanceC(int qpos, int vpos, double cap, double voltage) {
  double g, i;
  conductor(cap, g);
  addY(qpos, vpos, +g);
  i = pol * (g * voltage);
  addI(qpos, +i);
}

// The function initializes the histories of a circuit having the given age.
void circuit::initHistory(double age) {
  nHistories = getSize() + getVoltageSources();
  histories = new history[nHistories];
  setHistoryAge(age);
}

// Sets the age of all circuit histories
void circuit::setHistoryAge(double age) {
  for (int i = 0; i < nHistories; i++) {
    histories[i].setAge(age);
  }
}

// The function deletes the histories for the transient analysis.
void circuit::deleteHistory() {
  if (histories != nullptr) {
    delete[] histories;
    histories = nullptr;
  }
  setHistory(false);
}

// Truncates the transient analysis history (i.e. removes values newer
// newer than time tcut).
void circuit::truncateHistory(double tcut) {
  if (histories != nullptr) {
    for (int i = 0; i < nHistories; i++) {
      histories[i].truncate(tcut);
    }
  }
}

// Appends a history value.
void circuit::appendHistory(int n, double val) { histories[n].push_back(val); }

// Returns the required age of the history.
double circuit::getHistoryAge() {
  if (histories)
    return histories[0].getAge();
  return 0.0;
}

// Returns size of the history
int circuit::getHistorySize() { return histories[0].size(); }

// Returns the time with the specified index
double circuit::getHistoryTFromIndex(int idx) { return histories[0].getTfromidx(idx); }

/* This function should be used to apply the time vector history to
   the value histories of a circuit. */
void circuit::applyHistory(history *h) {
  for (int i = 0; i < nHistories; i++) {
    histories[i].apply(*h);
  }
}

// Returns voltage at the given time for the given node.
double circuit::getV(int port, double t) { return histories[port].nearest(t); }

// Returns voltage at the given index from the history for the given node.
double circuit::getV(int port, int idx) { return histories[port].getValfromidx(idx); }

// Returns current at the given time for the given voltage source.
double circuit::getJ(int nr, double t) { return histories[nr + getSize()].nearest(t); }

} // namespace qucs
