/*
 * mslange.h - parallel coupled microstrip lines class definitions
 *
 * Copyright (C) 2004, 2005, 2006, 2008 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __MSLANGE_H__
#define __MSLANGE_H__

class mslange : public qucs::circuit
{
 public:
  CREATOR (mslange);
  void initDC (void);
  void calcSP (double);
  void calcNoiseSP (double);
  void calcPropagation (double);
  void initAC (void);
  void calcAC (double);
  void calcNoiseAC (double);
  void saveCharacteristics (double);

  static void analysQuasiStatic (double, double, double,
				 double, double, const char * const,
				 double&, double&, double&,
				 double&);
  static void analyseDispersion (double, double, double,
				 double, double, double,
				 double, double, double, const char * const,
				 double&, double&, double&,
				 double&);

 private:
  double ae, be, ze, ao, bo, zo, ee, eo;
};

#endif /* __MSLANGE_H__ */
