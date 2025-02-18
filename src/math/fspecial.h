/*
 * fspecial.h - special functions definitions
 *
 * Copyright (C) 2006 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __FSPECIAL_H__
#define __FSPECIAL_H__

namespace fspecial {

//#ifndef HAVE_ERF
  double     erf (double);
//#endif
//#ifndef HAVE_ERFC
  double    erfc (double);
//#endif
  double  erfinv (double);
  double erfcinv (double);
  double ltqnorm (double);
  double      i0 (double);

  void        ellip_ke (double, double&, double&);
  double ellip_rf (double, double, double);
  double ellip_sncndn (double, double,
			    double&, double&, double&);

} // namespace

#endif /* __FSPECIAL_H__ */
