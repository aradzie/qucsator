/*
 * sweep.h - variable sweep class definitions
 *
 * Copyright (C) 2004, 2008 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __SWEEP_H__
#define __SWEEP_H__

#include <string>

namespace qucs {

enum sweep_type {
  SWEEP_UNKNOWN = -1,
  SWEEP_CONSTANT,
  SWEEP_LINEAR,
  SWEEP_LOGARITHMIC,
  SWEEP_LIST
};

class object;

class sweep : public object {
public:
  sweep();
  explicit sweep(const std::string &);
  sweep(const sweep &) = delete;
  ~sweep();
  int getSize() const { return this->size; }
  int getType() const { return this->type; }
  void setSize(int);
  double get(int);
  void set(int, double);
  double next();
  double prev();
  void reverse();
  void reset() { this->counter = 0; };
  object *getParent() const { return this->parent; }
  void setParent(object *parent) { this->parent = parent; }
  char *toString();

protected:
  int type;

private:
  double *data;
  int size;
  char *text;
  int counter;
  object *parent;
};

class linsweep : public sweep {
public:
  linsweep();
  explicit linsweep(const std::string &);
  ~linsweep();
  void create(double, double, int);
};

class logsweep : public sweep {
public:
  logsweep();
  explicit logsweep(const std::string &);
  ~logsweep();
  void create(double, double, int);
};

class consweep : public sweep {
public:
  consweep();
  explicit consweep(const std::string &);
  ~consweep();
  void create(double);
};

class lstsweep : public sweep {
public:
  lstsweep();
  explicit lstsweep(const std::string &);
  ~lstsweep();
  void create(int);
};

} // namespace qucs

#endif /* __SWEEP_H__ */
