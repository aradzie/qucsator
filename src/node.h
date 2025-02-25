/*
 * node.h - node class definitions
 *
 * Copyright (C) 2003, 2004 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __NODE_H__
#define __NODE_H__

#include <string>

namespace qucs {

class circuit;

class node final {
public:
  node() : name(), port(0), internal(0), cir(nullptr), nNode(0) {}
  node(char *const n) : name(n), port(0), internal(0), cir(nullptr), nNode(0) {}

  void setName(const std::string &n) { this->name = n; }
  const char *getName() const { return this->name.c_str(); }

  void setPort(const int p) { this->port = p; }
  int getPort() const { return this->port; }

  void setCircuit(circuit *const c) { this->cir = c; }
  circuit *getCircuit() const { return this->cir; }

  void setInternal(int i) { internal = i; }
  int getInternal() { return internal; }

  void setNode(const int n) { this->nNode = n; }
  int getNode() const { return this->nNode; }

private:
  std::string name; // Name, like "gnd", "_net2", etc.
  int port;         // The port number of this node.
  int internal;
  circuit *cir;
  int nNode; // The unique number of this node. Only used in the HB solver.
};

} // namespace qucs

#endif /* __NODE_H__ */
