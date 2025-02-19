/*
 * integrator.h - integrator class definitions
 *
 * Copyright (C) 2004, 2005, 2006 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

#include "states.h"

#define MODE_NONE 0
#define MODE_INIT 1

namespace qucs {

class integrator : public states<double> {
public:
  integrator();
  integrator(const integrator &);
  ~integrator();

  typedef void (*integrate_func_t)(integrator *, int, double, double &, double &);
  void setIntegration(integrate_func_t f) { integrate_func = f; }
  typedef void (*conductor_func_t)(integrator *, double, double &);
  void setConductance(conductor_func_t f) { conductor_func = f; }
  void integrate(int, double, double &, double &);
  void conductor(double, double &);
  void setOrder(int o) { order = o; }
  int getOrder() { return order; }
  void setMode(int s) { state = s; }
  int getMode() { return state; }
  void setCoefficients(double *c) { coefficients = c; }
  double *getCoefficients() { return coefficients; }

private:
  int order;
  int state;
  double *coefficients;
  integrate_func_t integrate_func;
  conductor_func_t conductor_func;
};

} // namespace qucs

#endif /* __INTEGRATOR_H__ */
