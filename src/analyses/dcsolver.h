/*
 * dcsolver.h - DC solver class definitions
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

#ifndef __DCSOLVER_H__
#define __DCSOLVER_H__

#include "nasolver.h"

namespace qucs {

class dcsolver final : public nasolver<double> {
public:
  ACREATOR(dcsolver);
  explicit dcsolver(char *);
  dcsolver(const dcsolver &) = delete;
  ~dcsolver() override;
  int solve() override;

private:
  void initDC();
  static void calcDC(dcsolver *);
  void saveOperatingPoints();

private:
  int saveOPs;
};

} // namespace qucs

#endif /* __DCSOLVER_H__ */
