/*
 * opamp.h - operational amplifier class definitions
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

#ifndef __OPAMP_H__
#define __OPAMP_H__

class opamp final : public qucs::circuit {
public:
  CREATOR(opamp);
  void calcDC();
  void calcTR(double);
  void initAC();
  void initDC();
  void initSP();
  void initTR();
  void calcOperatingPoints();

private:
  double gv;
};

#endif /* __OPAMP_H__ */
