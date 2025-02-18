/*
 * mosfet.h - mosfet class definitions
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

#ifndef __MOSFET_H__
#define __MOSFET_H__

class mosfet : public qucs::circuit
{
 public:
  CREATOR (mosfet);
  void calcSP (double);
  void calcNoiseSP (double);
  void calcDC (void);
  void initDC (void);
  void restartDC (void);
  void initModel (void);
  void saveOperatingPoints (void);
  void calcOperatingPoints (void);
  void loadOperatingPoints (void);
  void initAC (void);
  void calcAC (double);
  void calcNoiseAC (double);
  void initTR (void);
  void calcTR (double);

 private:
  double transientChargeTR (int, double&, double, double);
  double transientChargeSR (int, double&, double, double);
  qucs::matrix calcMatrixY (double);
  qucs::matrix calcMatrixCy (double);

 private:
  double UbsPrev, UbdPrev, UgsPrev, UgdPrev, UdsPrev, Udsat, Uon;
  double gbs, gbd, gm, gds, gmb, Ids, DrainControl, SourceControl;
  double Leff, MOSdir, beta, Cox, Phi, Ga, Vto, Rs, Rd;
  double Qgd, Qgs, Qbd, Qbs, Qgb, Ibs, Ibd;
  double Ugd, Ugs, Ubs, Ubd, Uds, Ugb;
  int transientMode;
  qucs::circuit * rs;
  qucs::circuit * rd;
  qucs::circuit * rg;
};

#endif /* __MOSFET_H__ */
