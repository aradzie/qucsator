/*
 * acsolver.h - AC solver class definitions
 *
 * Copyright (C) 2004, 2005, 2007, 2008 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __ACSOLVER_H__
#define __ACSOLVER_H__

#include "nasolver.h"

namespace qucs {

class sweep;
class vector;

class acsolver final : public nasolver<nr_complex_t> {
public:
  ACREATOR(acsolver);
  explicit acsolver(char *);
  acsolver(const acsolver &) = delete;
  ~acsolver() override;
  int solve() override;
  void solve_noise();
  static void calc(acsolver *);
  void init();
  void saveAllResults(double);
  void saveNoiseResults(qucs::vector *);

private:
  sweep *swp;
  double freq;
  int noise;
  tvector<double> *xn;
};

} // namespace qucs

#endif /* __ACSOLVER_H__ */
