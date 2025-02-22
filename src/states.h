/*
 * states.h - save-state variable class definitions
 *
 * Copyright (C) 2004 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __STATES_H__
#define __STATES_H__

namespace qucs {

/*
 * This class is used for storing sets of states for use
 * by the transient integrators.
 */
template <class state_type_t> class states {
public:
  states();
  states(const states &);
  ~states();

  state_type_t getState(int, int n = 0);
  void setState(int, state_type_t, int n = 0);
  void initStates();
  void clearStates();
  int getStates() { return nstates; }
  void setStates(int n) { nstates = n; }
  void nextState();
  void prevState();
  void fillState(int, state_type_t);
  void saveState(int, state_type_t *);
  void inputState(int, state_type_t *);

private:
  // Array for holding all the sets of states.
  // Multiple sets of states are stored in one large array which is indexed appropriately
  // to get the right state set and value.
  state_type_t *stateval;
  // The number of sets of states stored.
  int nstates;
  int currentstate;
};

} // namespace qucs

#include "states.cpp"

#endif /* __STATES_H__ */
