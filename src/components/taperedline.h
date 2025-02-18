/*
 * taperedline.h - ideal tapered transmission line class definition
 *
 * Copyright (C) 2015 Claudio Girardi <in3otd@qsl.net>
 * Copyright (C) 2015 Andres Martinez-Mera <andresmartinezmera@gmail.com>
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

#ifndef TAPEREDLINE_H
#define TAPEREDLINE_H
#include "matrix.h"

const int Nsteps = 20; // Number of sections used to approximate the taper

class taperedline : public qucs::circuit
{
 public:
  CREATOR (taperedline);
  void calcSP (double);
  void initDC (void);
  void initAC (void);
  void initSP (void);
  void calcAC (double);
  void calcNoiseAC (double);
  void calcNoiseSP (double);
private:
  void calcABCDparams(double);
  void calcImpedanceProfile();
  double calcExponential(double, double, double, double);
  double calcLinear(double, double, double, double);
  double calcTriangular(double, double, double, double);
  double calcKlopfenstein(double, double, double, double, double);
  double phi(double, double);
  qucs::matrix ABCD;
  double Zprofile[Nsteps];
};

#endif /* __taperedline_H__ */
