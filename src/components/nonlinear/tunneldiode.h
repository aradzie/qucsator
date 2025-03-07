/*
 * tunneldiode.h - resonance tunnel diode class definitions
 *
 * Copyright (C) 2011 Michael Margraf <michael.margraf@alumni.tu-berlin.de>
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

#ifndef __TUNNELDIODE_H__
#define __TUNNELDIODE_H__

class tunneldiode final : public qucs::circuit {
public:
  CREATOR(tunneldiode);
  void calcAC(double) override;
  void calcDC() override;
  void calcSP(double) override;
  void calcTR(double) override;
  void initAC() override;
  void initDC() override;
  void initTR() override;
  void saveOperatingPoints() override;
  void loadOperatingPoints();
  void calcOperatingPoints() override;

private:
  double Ud, gd, Id, Qd;

private:
  qucs::matrix calcMatrixY(double);

  void calcId(double, double &, double &);
};

#endif /* __TUNNELDIODE_H__ */
