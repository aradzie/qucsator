/*
 * analysis.h - analysis class definitions
 *
 * Copyright (C) 2003, 2004, 2005, 2006, 2007 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include "complex.h"
#include "object.h"
#include "ptrlist.h"

#define SAVE_OPS 1 // save operating points
#define SAVE_ALL 2 // also save sub-circuit nodes and operating points
#define SAVE_CVS 4 // save characteristic values

#define ACREATOR(val)                                                                              \
  val();                                                                                           \
  static analysis *create() { return new val(); }                                                  \
  static struct define_t anadef;                                                                   \
  static struct define_t *definition() { return &anadef; }

namespace qucs {

class dataset;
class net;
class environment;
class sweep;
class vector;

enum analysis_type {
  ANALYSIS_UNKNOWN = -1,
  ANALYSIS_SWEEP,
  ANALYSIS_DC,
  ANALYSIS_AC,
  ANALYSIS_HBALANCE,
  ANALYSIS_TRANSIENT,
  ANALYSIS_SPARAMETER,
  ANALYSIS_E_TRANSIENT
};

class analysis : public object {
public:
  analysis();
  explicit analysis(const std::string &);
  analysis(const analysis &) = delete;
  virtual ~analysis();
  virtual int solve() { return 0; }
  virtual int initialize() { return 0; }
  virtual int cleanup() { return 0; }
  virtual bool isExternal() { return false; }
  dataset *getData() const { return this->data; }
  void setData(dataset *data) { this->data = data; }
  net *getNet() { return this->subnet; }
  void setNet(net *netlist) { this->subnet = netlist; }
  environment *getEnv() { return this->env; }
  void setEnv(environment *env) { this->env = env; }
  ptrlist<analysis> *getAnalysis() { return this->actions; }
  void setAnalysis(ptrlist<analysis> *actions) { this->actions = actions; }
  void addAnalysis(analysis *);
  void delAnalysis(analysis *);
  int getType() { return this->type; }
  void setType(int type) { this->type = type; }
  sweep *createSweep(const std::string &);
  void saveVariable(const std::string &, nr_complex_t, qucs::vector *);
  bool getProgress() const { return this->progress; }
  void setProgress(bool progress) { this->progress = progress; }

protected:
  int runs;
  int type;
  net *subnet;
  dataset *data;
  environment *env;
  ptrlist<analysis> *actions;
  bool progress;
};

} // namespace qucs

#endif /* __ANALYSIS_H__ */
