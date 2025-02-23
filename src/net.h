/*
 * net.h - net class definitions
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

#ifndef __NET_H__
#define __NET_H__

#include <string>

#include "ptrlist.h"

namespace qucs {

class circuit;
class node;
class nodelist;
class nodeset;
class analysis;
class dataset;
class environment;

class net : public object {
public:
  explicit net();
  explicit net(const std::string &);
  net(const net &) = delete;
  ~net();
  circuit *getRoot() { return root; }
  void setRoot(circuit *c) { root = c; }
  void insertCircuit(circuit *);
  void removeCircuit(circuit *, int dropping = 1);
  int containsCircuit(circuit *);
  int checkCircuitChain();
  void list();
  void reducedCircuit(circuit *);
  node *findConnectedNode(node *);
  node *findConnectedCircuitNode(node *);
  void insertedCircuit(circuit *);
  void insertedNode(node *);
  void insertAnalysis(analysis *);
  void removeAnalysis(analysis *);
  dataset *runAnalysis(int &);
  void getDroppedCircuits(nodelist *nodes = nullptr);
  void deleteUnusedCircuits(nodelist *nodes = nullptr);
  int getPorts() { return nPorts; }
  int getReduced() { return reduced; }
  void setReduced(int r) { reduced = r; }
  int getVoltageSources() { return nSources; }
  void setVoltageSources(int n) { nSources = n; }
  analysis *findAnalysis(const std::string &) const;
  analysis *findAnalysis(int);
  analysis *findSecondOrder();
  analysis *getChildAnalysis(analysis *);
  const char *getChild(analysis *) const;
  void orderAnalysis();
  analysis *findLastOrder(analysis *);
  ptrlist<analysis> *findLastOrderChildren(analysis *);
  void sortChildAnalyses(analysis *);
  int containsAnalysis(analysis *, int);
  environment *getEnv() { return env; }
  void setEnv(environment *e) { env = e; }
  int countPorts();
  int countNodes();
  int isNonLinear();
  void addNodeset(nodeset *);
  void delNodeset();
  nodeset *getNodeset() { return nset; }
  void setSrcFactor(double f) { srcFactor = f; }
  double getSrcFactor() { return srcFactor; }
  void setActionNetAll(net *);

private:
  nodeset *nset;
  circuit *drop;
  circuit *root;
  ptrlist<analysis> *actions;
  ptrlist<analysis> *orgacts;
  environment *env;
  int nPorts;
  int nSources; // ARA: The number of voltage sources in the whole circuit. This influences the size of the system of linear equations.
  int nCircuits;
  int reduced;
  int inserted;
  int insertedNodes;
  double srcFactor;
};

} // namespace qucs

#endif /* __NET_H__ */
