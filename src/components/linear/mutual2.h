/*
 * mutual2.h - three mutual inductors class definitions
 *
 * Copyright (C) 2005, 2008 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __MUTUAL2_H__
#define __MUTUAL2_H__

class mutual2 final : public qucs::circuit {
public:
  CREATOR(mutual2);
  void calcSP(double) override;
  void initDC() override;
  void initAC() override;
  void calcAC(double) override;
  void initTR() override;
  void calcTR(double) override;

private:
  qucs::matrix calcMatrixY(double);
};

#endif /* __MUTUAL2_H__ */
