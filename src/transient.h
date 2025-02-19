/*
 * transient.h - transient helper class definitions
 *
 * Copyright (C) 2004, 2006 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __TRANSIENT_H__
#define __TRANSIENT_H__

namespace qucs {

enum integrator_type {
  INTEGRATOR_UNKNOWN = -1,
  INTEGRATOR_EULER = 0,
  INTEGRATOR_TRAPEZOIDAL = 1,
  INTEGRATOR_GEAR = 2,
  INTEGRATOR_ADAMSMOULTON = 3,
  INTEGRATOR_ADAMSBASHFORD = 4
};

class circuit;
class integrator;

namespace transient {

void calcCorrectorCoeff(int, int, double *, double *);
void calcPredictorCoeff(int, int, double *, double *);
void getConductance(integrator *, double, double &);
void integrateEuler(integrator *, int, double, double &, double &);
void integrateBilinear(integrator *, int, double, double &, double &);
void integrateGear(integrator *, int, double, double &, double &);
void integrateMoulton(integrator *, int, double, double &, double &);
void setIntegrationMethod(circuit *, int);
int correctorType(const char *const, int &);
int correctorType(int, int);
int predictorType(int, int, int &);
double getCorrectorError(int, int);
double getPredictorError(int, int);

} // namespace transient

} // namespace qucs

#endif /* __TRANSIENT_H__ */
