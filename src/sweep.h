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
  sweep(const std::string &);
  sweep(sweep &);
  ~sweep();
  int getSize() { return size; }
  int getType() { return type; }
  double get(int);
  double next();
  double prev();
  void set(int, double);
  void setSize(int);
  char *toString();
  void reverse();
  void reset() { counter = 0; };
  object *getParent() { return parent; }
  void setParent(object *p) { parent = p; }

protected:
  int type;

private:
  double *data;
  int size;
  char *txt;
  int counter;
  object *parent;
};

class linsweep : public sweep {
public:
  linsweep();
  linsweep(const std::string &);
  ~linsweep();
  void create(double, double, int);
};

class logsweep : public sweep {
public:
  logsweep();
  logsweep(const std::string &);
  ~logsweep();
  void create(double, double, int);
};

class consweep : public sweep {
public:
  consweep();
  consweep(const std::string &);
  ~consweep();
  void create(double);
};

class lstsweep : public sweep {
public:
  lstsweep();
  lstsweep(const std::string &);
  ~lstsweep();
  void create(int);
};

} // namespace qucs

#endif /* __SWEEP_H__ */
