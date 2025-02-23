/*
 * circuit.h - circuit class definitions
 *
 * Copyright (C) 2003-2008 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __CIRCUIT_H__
#define __CIRCUIT_H__

#include "complex.h"
#include "object.h"
#include "pair.h"

#define NODE_1 0
#define NODE_2 1
#define NODE_3 2
#define NODE_4 3
#define NODE_5 4
#define NODE_6 5
#define VSRC_1 0
#define VSRC_2 1
#define VSRC_3 2
#define VSRC_4 3
#define VSRC_5 4

#define MODFLAG(val, bit)                                                                          \
  if (val)                                                                                         \
    flag |= (bit);                                                                                 \
  else                                                                                             \
    flag &= ~(bit);
#define RETFLAG(bit) ((flag & (bit)) != 0)

#define CREATOR(val)                                                                               \
  val();                                                                                           \
  static qucs::circuit *create() { return new val(); }                                             \
  static struct define_t cirdef;                                                                   \
  static struct define_t *definition() { return &cirdef; }

#include <map>
#include <string>

#include "integrator.h"
#include "valuelist.h"

namespace qucs {

enum circuit_flag {
  CIRCUIT_ENABLED = 1,
  CIRCUIT_LINEAR = 2,
  CIRCUIT_ORIGINAL = 4,
  CIRCUIT_VSOURCE = 8,
  CIRCUIT_ISOURCE = 16,
  CIRCUIT_INTVSOURCE = 32,
  CIRCUIT_VARSIZE = 64,
  CIRCUIT_PROBE = 128,
  CIRCUIT_HISTORY = 256,
};

class node;
class property;
class substrate;
class matrix;
class net;
class environment;
class history;

/* The base class for all circuit elements and provides the functionality
 * required for all simulation types. */
class circuit : public object, public integrator {
public:
  circuit *getNext() const { return this->next; }
  void setNext(circuit *const o) { this->next = o; }
  circuit *getPrev() const { return prev; }
  void setPrev(circuit *const o) { this->prev = o; }

private:
  circuit *next;
  circuit *prev;

public:
  explicit circuit();
  explicit circuit(int);
  circuit(const circuit &) = delete;
  virtual ~circuit();

  virtual void initDC() { allocMatrixMNA(); }
  virtual void calcDC() {}
  virtual void restartDC() {} // ARA: What is this method for?

  virtual void initTR() { allocMatrixMNA(); }
  virtual void calcTR(double) {}

  virtual void initAC() { allocMatrixMNA(); }
  virtual void calcAC(double) {}
  virtual void initNoiseAC() { allocMatrixN(vsources); }
  virtual void calcNoiseAC(double) {}

  virtual void initSP() { allocMatrixS(); }
  virtual void calcSP(double) {}
  virtual void initNoiseSP() { allocMatrixN(); }
  virtual void calcNoiseSP(double) {}

  virtual void initHB() { allocMatrixMNA(); }
  virtual void initHB(int) { allocMatrixMNA(); }
  virtual void calcHB(double) {}
  virtual void calcHB(int) {}

  virtual void calcOperatingPoints() {}
  virtual void saveOperatingPoints() {}

  virtual void calcCharacteristics(double) {}
  virtual void saveCharacteristics(double) {}
  virtual void saveCharacteristics(nr_complex_t) {}

  void setNode(int, const std::string &, int intern = 0);
  node *getNode(int) const;

  int getType() const { return type; }

  /* Gets the number of ports the circuit element has. */
  int getSize() const { return size; }
  /* Sets/changes the number of ports the circuit element has.
   * On setting this value, previously stored node and matrix
   * information is completely lost unless the new size equals
   * the original size. */
  void setSize(int);

  void setEnabled(bool e) { MODFLAG(e, CIRCUIT_ENABLED); }
  bool isEnabled() const { return RETFLAG(CIRCUIT_ENABLED); }

  void setVariableSized(bool v) { MODFLAG(v, CIRCUIT_VARSIZE); }
  bool isVariableSized() const { return RETFLAG(CIRCUIT_VARSIZE); }

  void setProbe(bool p) { MODFLAG(p, CIRCUIT_PROBE); }
  bool isProbe() const { return RETFLAG(CIRCUIT_PROBE); }

  void setNet(net *n) { subnet = n; }
  net *getNet() const { return subnet; }

  // subcircuitry

  void setSubcircuit(const std::string &n) { subcircuit = n; };
  std::string getSubcircuit() const { return subcircuit; }

  // environment specific

  void setEnv(environment *e) { env = e; }
  environment *getEnv() const { return env; }

  // nodal analyses helpers

  void setInternalVoltageSource(bool i) { MODFLAG(i, CIRCUIT_INTVSOURCE); }
  bool isInternalVoltageSource() const { return RETFLAG(CIRCUIT_INTVSOURCE); }

  void setVoltageSource(int s) { vsource = s; }
  int getVoltageSource() const { return vsource; }

  void setVoltageSources(int s) { vsources = s; }
  int getVoltageSources() const { return vsources; }

  void voltageSource(int, int, int, double value = 0.0);

  void setVSource(bool v) { MODFLAG(v, CIRCUIT_VSOURCE); }
  bool isVSource() const { return RETFLAG(CIRCUIT_VSOURCE); }

  void setISource(bool i) { MODFLAG(i, CIRCUIT_ISOURCE); }
  bool isISource() const { return RETFLAG(CIRCUIT_ISOURCE); }

  /* Sets the number of internal noise voltage sources for AC analysis. */
  void setNoiseSources(int s) { nsources = s; }
  /* Returns the number of internal noise sources for AC analysis. */
  int getNoiseSources() const { return nsources; }

  // transient analyses helpers

  void transientCapacitance(int, int, int, double, double, double);
  void transientCapacitance(int, int, double, double, double);
  void transientCapacitanceQ(int, int, int, double);
  void transientCapacitanceQ(int, int, double);
  void transientCapacitanceC(int, int, int, int, double, double);
  void transientCapacitanceC(int, int, double, double);
  void transientCapacitanceC2V(int, int, int, double, double);
  void transientCapacitanceC2Q(int, int, int, double, double);

  void setDelta(double *d) { deltas = d; }
  double *getDelta() const { return deltas; }

  // history specific functionality

  bool hasHistory() const { return RETFLAG(CIRCUIT_HISTORY); }
  void setHistory(bool h) { MODFLAG(h, CIRCUIT_HISTORY); }
  void initHistory(double);
  void deleteHistory();
  void truncateHistory(double);
  void appendHistory(int, double);
  void applyHistory(history *);

  double getV(int, double);
  double getV(int, int);
  double getJ(int, double);

  double getHistoryAge();
  void setHistoryAge(double);
  int getHistorySize();
  double getHistoryTFromIndex(int);

  // s-parameter helpers

  int getPort() const { return pacport; }
  void setPort(int p) { pacport = p; }
  int getInserted() const { return inserted; }
  void setInserted(int i) { inserted = i; }
  bool isOriginal() const { return RETFLAG(CIRCUIT_ORIGINAL); }
  void setOriginal(bool o) { MODFLAG(o, CIRCUIT_ORIGINAL); }

  // microstrip helpers

  /* Gets the current substrate of the circuit object. Used for microstrip components only. */
  substrate *getSubstrate() const { return subst; };
  /* Sets the substrate of the circuit object. */
  void setSubstrate(substrate *s) { subst = s; };

  // matrix entry modificators

  /* Returns the S-parameter at the given matrix position. */
  nr_complex_t getS(int x, int y) const { return MatrixS[y + x * size]; }
  /* Returns the noise-correlation-parameter at the given matrix position. */
  nr_complex_t getN(int r, int c) const { return MatrixN[c + r * (size + nsources)]; };
  /* Returns the circuits G-MNA matrix value depending on the port numbers. */
  nr_complex_t getY(int r, int c) const { return MatrixY[r * size + c]; }
  /* Returns the circuits B-MNA matrix value of the given voltage source built in the circuit
   * depending on the port number. */
  nr_complex_t getB(int port, int nr) const { return MatrixB[(nr - vsource) * size + port]; }
  /* Returns the circuits C-MNA matrix value of the given voltage source built in the circuit
   * depending on the port number. */
  nr_complex_t getC(int nr, int port) const { return MatrixC[(nr - vsource) * size + port]; }
  /* Returns the circuits D-MNA matrix value of the given voltage source built in the circuit. */
  nr_complex_t getD(int r, int c) const { return MatrixD[(r - vsource) * vsources + c - vsource]; }
  /* Returns the circuits C-HB matrix value depending on the port numbers. */
  nr_complex_t getQV(int r, int c) const { return MatrixQV[r * size + c]; }
  /* Returns the circuits GV-HB vector value depending on the port number. */
  nr_complex_t getGV(int port) const { return VectorGV[port]; }
  /* Returns the circuits CV-HB vector value depending on the port number. */
  nr_complex_t getCV(int port) const { return VectorCV[port]; }
  /* Returns the circuits E-MNA matrix value of the given voltage source built in the circuit. */
  nr_complex_t getE(int nr) const { return VectorE[nr - vsource]; }
  /* Returns the circuits I-MNA matrix value of the current source built in the circuit. */
  nr_complex_t getI(int port) const { return VectorI[port]; }
  /* Returns the circuits J-MNA matrix value of the given voltage source built in the circuit. */
  nr_complex_t getJ(int nr) const { return VectorJ[nr]; }
  /* Returns the circuits voltage value at the given port. */
  nr_complex_t getV(int port) const { return VectorV[port]; }
  /* Returns the circuits Q-HB vector value. */
  nr_complex_t getQ(int port) const { return VectorQ[port]; }

  /* Sets the S-parameter at the given matrix position. */
  void setS(int x, int y, nr_complex_t z) { MatrixS[y + x * size] = z; }
  /* Sets the noise-correlation-parameter at the given matrix position. */
  void setN(int r, int c, nr_complex_t z) { MatrixN[c + r * (size + nsources)] = z; }
  /* Sets the circuits G-MNA matrix value depending on the port numbers. */
  void setY(int r, int c, nr_complex_t y) { MatrixY[r * size + c] = y; }
  /* Modifies the circuits G-MNA matrix value depending on the port numbers. */
  void addY(int r, int c, nr_complex_t y) { MatrixY[r * size + c] += y; }
  /* Modifies the circuits G-MNA matrix value depending on the port numbers. */
  void addY(int r, int c, double y) { MatrixY[r * size + c] += y; }
  /* Sets the circuits B-MNA matrix value of the given voltage source built in the circuit depending
   * on the port number. */
  void setB(int port, int nr, nr_complex_t z) { MatrixB[nr * size + port] = z; }
  /* Sets the circuits C-MNA matrix value of the given voltage source built in the circuit depending
   * on the port number. */
  void setC(int nr, int port, nr_complex_t z) { MatrixC[nr * size + port] = z; }
  /* Sets the circuits D-MNA matrix value of the given voltage source built in the circuit. */
  void setD(int r, int c, nr_complex_t z) { MatrixD[r * vsources + c] = z; }
  /* Sets the circuits C-HB matrix value depending on the port numbers. */
  void setQV(int r, int c, nr_complex_t qv) { MatrixQV[r * size + c] = qv; }
  /* Sets the circuits GV-HB matrix value depending on the port number. */
  void setGV(int port, nr_complex_t gv) { VectorGV[port] = gv; }
  /* Sets the circuits CV-HB matrix value depending on the port number. */
  void setCV(int port, nr_complex_t cv) { VectorCV[port] = cv; }
  /* Sets the circuits E-MNA matrix value of the given voltage source built in the circuit. */
  void setE(int nr, nr_complex_t z) { VectorE[nr] = z; }
  /* Sets the circuits I-MNA matrix value of the current source built in the circuit depending on
   * the port number. */
  void setI(int port, nr_complex_t z) { VectorI[port] = z; }
  /* Modifies the circuits I-MNA matrix value of the current source built in the circuit depending
   * on the port number. */
  void addI(int port, nr_complex_t i) { VectorI[port] += i; }
  /* Same as above with different argument type. */
  void addI(int port, double i) { VectorI[port] += i; }
  /* Sets the circuits J-MNA matrix value of the given voltage source built in the circuit. */
  void setJ(int nr, nr_complex_t z) { VectorJ[nr - vsource] = z; }
  /* Sets the circuits voltage value at the given port. */
  void setV(int port, nr_complex_t z) { VectorV[port] = z; }
  /* Sets the circuits Q-HB vector value. */
  void setQ(int port, nr_complex_t q) { VectorQ[port] = q; }
  /* Sets the circuits G-MNA matrix value depending on the port numbers. */
  void setG(int r, int c, double y) { MatrixY[r * size + c] = y; }

  void clearB() { memset(MatrixB, 0, sizeof(nr_complex_t) * size * vsources); }
  void clearC() { memset(MatrixC, 0, sizeof(nr_complex_t) * size * vsources); }
  void clearD() { memset(MatrixD, 0, sizeof(nr_complex_t) * vsources * vsources); }
  void clearE() { memset(VectorE, 0, sizeof(nr_complex_t) * vsources); }
  void clearJ() { memset(VectorJ, 0, sizeof(nr_complex_t) * vsources); }
  void clearI() { memset(VectorI, 0, sizeof(nr_complex_t) * size); }
  void clearV() { memset(VectorV, 0, sizeof(nr_complex_t) * size); }
  void clearY() { memset(MatrixY, 0, sizeof(nr_complex_t) * size * size); }

  // operating point functionality

  void addOperatingPoint(const std::string &name, double);
  double getOperatingPoint(const std::string &name);
  void setOperatingPoint(const std::string &name, double);
  int hasOperatingPoint(const std::string &name);
  valuelist<pair> &getOperatingPoints() { return oper; }

  // characteristics functionality

  void addCharacteristic(const std::string &name, double);
  double getCharacteristic(const std::string &name);
  void setCharacteristic(const std::string &name, double);
  int hasCharacteristic(const std::string &name);
  valuelist<pair> &getCharacteristics() { return charac; }

  // differentiate between linear and non-linear circuits

  void setNonLinear(bool l) { MODFLAG(!l, CIRCUIT_LINEAR); }
  bool isNonLinear() const { return !RETFLAG(CIRCUIT_LINEAR); }

  // miscellaneous functionality

  static std::string createInternal(const std::string &, const std::string &);
  void setInternalNode(int, const std::string &);
  void print();

  // matrix operations

  void allocMatrixS();
  void allocMatrixN(int sources = 0);
  void allocMatrixMNA();
  void freeMatrixMNA();
  void allocMatrixHB();
  void freeMatrixHB();
  void setMatrixS(matrix);
  matrix getMatrixS();
  void setMatrixN(matrix);
  matrix getMatrixN();
  void setMatrixY(matrix);
  matrix getMatrixY();

  static const double z0;

protected:
  int type;
  int pol;

private:
  int size; // The number of ports/nodes.
  int pacport; // ARA: Something for S-parameters.
  int vsource; // ARA: Index of the first voltage source in the MNA matrix.
  int vsources; // ARA: The number of internal voltage sources.
  int nsources; // ARA: The number of internal voltage noise sources for AC analysis.
  int inserted;
  int flag;
  nr_complex_t *MatrixS;
  nr_complex_t *MatrixN;
  nr_complex_t *MatrixY;  // ARA: Devices update this vector. The `nasolver` class reads this value.
                          // AC analysis.
  nr_complex_t *MatrixB;  // ARA: Devices update this vector. The `nasolver` class reads this value.
                          // DC/TR analysis.
  nr_complex_t *MatrixC;  // ARA: Devices update this vector. The `nasolver` class reads this value.
                          // DC/TR analysis.
  nr_complex_t *MatrixD;  // ARA: Devices update this vector. The `nasolver` class reads this value.
                          // DC/TR analysis.
  nr_complex_t *VectorE;  // ARA: Device voltages. Devices update this vector. The `nasolver` class
                          // reads this value. DC/TR analysis.
  nr_complex_t *VectorI;  // ARA: Device currents. Devices update this vector. The `nasolver` class
                          // reads this value. DC/TR analysis.
  nr_complex_t *VectorV;  // ARA: The `nasolver` class updates this vector with the solved node
                          // voltages. Devices read this value. DC/TR analysis.
  nr_complex_t *VectorJ;  // ARA: The `nasolver` class updates this vector with the solved branch
                          // currents. Devices read this value. DC/TR analysis.
  nr_complex_t *VectorQ;  // ARA: For HB analysis only.
  nr_complex_t *MatrixQV; // ARA: For HB analysis only.
  nr_complex_t *VectorGV; // ARA: For HB analysis only.
  nr_complex_t *VectorCV; // ARA: For HB analysis only.
  std::string subcircuit;
  node *nodes;
  substrate *subst;
  valuelist<pair> oper;
  valuelist<pair> charac;
  net *subnet;
  environment *env;
  double *deltas;
  int nHistories;
  history *histories;
};

} // namespace qucs

// typedef to make it easier to set up our factory
typedef qucs::circuit *maker_t();
// function typdefs to make it easier to set up our factories
typedef qucs::circuit *creator_t();
typedef struct define_t *defs_t();

// our global factories defined in module.cpp
extern "C" {
extern std::map<std::string, creator_t *> factorycreate;
extern std::map<std::string, defs_t *> factorydef;
}

#endif /* __CIRCUIT_H__ */
