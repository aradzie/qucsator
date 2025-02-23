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
  node *getNode(int);

  int getType() { return type; }

  /* Gets the number of ports the circuit element has. */
  int getSize() { return size; }
  /* Sets/changes the number of ports the circuit element has.
   * On setting this value, previously stored node and matrix
   * information is completely lost unless the new size equals
   * the original size. */
  void setSize(int);

  /* Returns true if the circuit element is enabled or false
   * otherwise. */
  bool isEnabled() { return RETFLAG(CIRCUIT_ENABLED); }
  /* Sets the circuit element to be enabled or disabled. */
  void setEnabled(bool e) { MODFLAG(e, CIRCUIT_ENABLED); }

  bool isVariableSized() { return RETFLAG(CIRCUIT_VARSIZE); }
  void setVariableSized(bool v) { MODFLAG(v, CIRCUIT_VARSIZE); }

  bool isProbe() { return RETFLAG(CIRCUIT_PROBE); }
  void setProbe(bool p) { MODFLAG(p, CIRCUIT_PROBE); }

  void setNet(net *n) { subnet = n; }
  net *getNet() { return subnet; }

  // subcircuitry
  std::string getSubcircuit() { return subcircuit; }
  void setSubcircuit(const std::string &);

  // environment specific
  environment *getEnv() { return env; }
  void setEnv(environment *e) { env = e; }

  // nodal analyses helpers
  void setInternalVoltageSource(bool i) { MODFLAG(i, CIRCUIT_INTVSOURCE); }
  bool isInternalVoltageSource() { return RETFLAG(CIRCUIT_INTVSOURCE); }
  void setVoltageSource(int s) { vsource = s; }
  int getVoltageSource() { return vsource; }
  int getVoltageSources();
  void setVoltageSources(int);
  void voltageSource(int, int, int, double value = 0.0);

  bool isVSource() { return RETFLAG(CIRCUIT_VSOURCE); }
  void setVSource(bool v) { MODFLAG(v, CIRCUIT_VSOURCE); }

  bool isISource() { return RETFLAG(CIRCUIT_ISOURCE); }
  void setISource(bool i) { MODFLAG(i, CIRCUIT_ISOURCE); }

  int getNoiseSources();
  void setNoiseSources(int);

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
  double *getDelta() { return deltas; }

  // history specific functionality
  bool hasHistory() { return RETFLAG(CIRCUIT_HISTORY); }
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
  int getPort() { return pacport; }
  void setPort(int p) { pacport = p; }
  int getInserted() { return inserted; }
  void setInserted(int i) { inserted = i; }
  bool isOriginal() { return RETFLAG(CIRCUIT_ORIGINAL); }
  void setOriginal(bool o) { MODFLAG(o, CIRCUIT_ORIGINAL); }

  // microstrip helpers
  substrate *getSubstrate();
  void setSubstrate(substrate *);

  // matrix entry modificators
  nr_complex_t getS(int, int);
  nr_complex_t getN(int, int);
  nr_complex_t getY(int, int);
  nr_complex_t getB(int, int);
  nr_complex_t getC(int, int);
  nr_complex_t getD(int, int);
  nr_complex_t getQV(int, int);
  nr_complex_t getGV(int);
  nr_complex_t getCV(int);
  nr_complex_t getE(int);
  nr_complex_t getI(int);
  nr_complex_t getJ(int);
  nr_complex_t getV(int);
  nr_complex_t getQ(int);
  double getG(int, int);
  void setS(int, int, nr_complex_t);
  void setN(int, int, nr_complex_t);
  void setY(int, int, nr_complex_t);
  void setB(int, int, nr_complex_t);
  void setC(int, int, nr_complex_t);
  void setD(int, int, nr_complex_t);
  void setQV(int, int, nr_complex_t);
  void setGV(int, nr_complex_t);
  void setCV(int, nr_complex_t);
  void setE(int, nr_complex_t);
  void setI(int, nr_complex_t);
  void setJ(int, nr_complex_t);
  void setV(int, nr_complex_t);
  void setQ(int, nr_complex_t);
  void setG(int, int, double);
  void clearB();
  void clearC();
  void clearD();
  void clearE();
  void clearI();
  void clearJ();
  void clearV();
  void clearY();
  void addY(int, int, nr_complex_t);
  void addY(int, int, double);
  void addI(int, nr_complex_t);
  void addI(int, double);

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
  bool isNonLinear() { return !RETFLAG(CIRCUIT_LINEAR); }

  // miscellaneous functionality
  void print();
  static std::string createInternal(const std::string &, const std::string &);
  void setInternalNode(int, const std::string &);

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
  int size;
  int pacport;
  int vsource; // ARA: Index of the voltage source in the MNA matrix.
  int vsources;
  int nsources;
  int inserted;
  int flag;
  nr_complex_t *MatrixS;
  nr_complex_t *MatrixN;
  nr_complex_t *MatrixY; // ARA: Devices update this vector. The `nasolver` class reads this value. AC analysis.
  nr_complex_t *MatrixB; // ARA: Devices update this vector. The `nasolver` class reads this value. DC/TR analysis.
  nr_complex_t *MatrixC; // ARA: Devices update this vector. The `nasolver` class reads this value. DC/TR analysis.
  nr_complex_t *MatrixD; // ARA: Devices update this vector. The `nasolver` class reads this value. DC/TR analysis.
  nr_complex_t *VectorE; // ARA: Device currents. Devices update this vector. The `nasolver` class reads this value. DC/TR analysis.
  nr_complex_t *VectorI; // ARA: Device voltages. Devices update this vector. The `nasolver` class reads this value. DC/TR analysis.
  nr_complex_t *VectorV; // ARA: The `nasolver` class updates this vector with the solved node voltages. Devices read this value. DC/TR analysis.
  nr_complex_t *VectorJ; // ARA: The `nasolver` class updates this vector with the solved branch currents. Devices read this value. DC/TR analysis.
  nr_complex_t *VectorQ; // ARA: For HB analysis only.
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
