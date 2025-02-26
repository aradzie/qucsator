/*
 * devstates.h - device state class definitions
 *
 * Copyright (C) 2006 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __DEVSTATES_H__
#define __DEVSTATES_H__

namespace qucs {

class devstates {
public:
  devstates();
  devstates(int, int);
  ~devstates();

  int deviceState();
  void deviceState(int);
  int deviceStates();
  void deviceStates(int, int);
  double operator()(int) const;
  double &operator()(int);
  double deviceVar(int) const;
  double &deviceVar(int);

private:
  int nstates;
  int nvars;
  int nstate;
  double *states;
  double *pstate;
};

} // namespace qucs

#endif /* __DEVSTATES_H__ */
