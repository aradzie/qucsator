/*
 * digital.h - digital base class definitions
 *
 * Copyright (C) 2005, 2006 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __DIGITAL_H__
#define __DIGITAL_H__

class digital : public qucs::circuit {
public:
  digital();
  ~digital();
  void calcAC(double) override;
  void calcDC() override;
  void calcSP(double) override;
  void calcTR(double) override;
  void initAC() override;
  void initDC() override;
  void initSP() override;
  void initTR() override;
  void calcOperatingPoints() override;

protected:
  virtual void calcOutput() {}
  virtual void calcDerivatives() {}
  double getVin(int);
  double calcTransfer(int);
  double calcTransferX(int);
  double calcDerivative(int);
  double calcDerivativeX(int);

protected:
  double *g;
  double Vout, Veq, Tdelay;
  int i;
  bool delay;

private:
  void initDigital();
  void freeDigital();
};

#endif /* __DIGITAL_H__ */
