/*
 * bjt.h - bipolar junction transistor class definitions
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

#ifndef __BJT_H__
#define __BJT_H__

class bjt : public qucs::circuit
{
 public:
  CREATOR (bjt);
  void calcSP (double);
  void calcNoiseSP (double);
  void initSP (void);
  void calcDC (void);
  void initDC (void);
  void restartDC (void);
  void saveOperatingPoints (void);
  void loadOperatingPoints (void);
  void calcOperatingPoints (void);
  void initAC (void);
  void calcAC (double);
  void calcNoiseAC (double);
  void initTR (void);
  void calcTR (double);

 private:
  void initModel (void);
  void processCbcx (void);
  qucs::matrix calcMatrixY (double);
  qucs::matrix calcMatrixCy (double);
  void excessPhase (int, double&, double&);

 private:
  double Ucs, Ubx, Ube, Ubc, Uce, UbePrev, UbcPrev;
  qucs::circuit * re;
  qucs::circuit * rc;
  qucs::circuit * rb;
  qucs::circuit * cbcx;
  double dQbedUbc, dQbdUbe, dQbdUbc, If, Qb, Ir, It;
  double gbei, gben, gbci, gbcn, gitf, gitr, gif, gir, Rbb, Ibe;
  double Qbe, Qbci, Qbcx, Qcs;
  bool doTR;
};

#endif /* __BJT_H__ */
