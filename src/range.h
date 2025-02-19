/*
 * range.h - range class definitions
 *
 * Copyright (C) 2006, 2008 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __RANGE_H__
#define __RANGE_H__

namespace qucs {

class range {
public:
  range();
  range(double, double);
  range(char, double, double, char);
  range(const range &);
  ~range();
  bool inside(double);
  bool outside(double);
  double lo(void) { return l; }
  double hi(void) { return h; }
  char *toString(void);

private:
  char il;  // interval boundary
  double l; // lower bound of the value
  double h; // upper bound of the value
  char ih;  // interval boundary
  char *txt;
};

} // namespace qucs

#endif /* __RANGE_H__ */
