/*
 * inductor.h - inductor class definitions
 *
 * Copyright (C) 2003, 2004, 2006, 2008 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __INDUCTOR_H__
#define __INDUCTOR_H__

class inductor final : public qucs::circuit {
public:
  CREATOR(inductor);
  void calcSP(double) override;
  void initDC() override;
  void calcDC() override;
  void initAC() override;
  void calcAC(double) override;
  void initTR() override;
  void calcTR(double) override;
  void initHB() override;
  void calcHB(double) override;
};

#endif /* __INDUCTOR_H__ */
