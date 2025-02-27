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

/**
 * Acts like N circular buffers of values.
 * Each such buffer has a fixed length of 8 elements.
 * Each such buffer has a unique numeric id, starting from zero.
 */
class states {
public:
  states();
  states(const states &) = delete;
  ~states();

  void setStates(int n);
  [[nodiscard]] int getStates() const;

  void fillState(int, double);
  void setState(int, double, int offset = 0);
  [[nodiscard]] double getState(int, int offset = 0) const;

  void nextState();

private:
  int size;
  double *values;
  int pos;
};

} // namespace qucs

#endif /* __STATES_H__ */
