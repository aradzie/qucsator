/*
 * device.h - device class definitions
 *
 * Copyright (C) 2004, 2005, 2006 Stefan Jahn <stefan@lkcc.org>
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

#ifndef __DEVICE_H__
#define __DEVICE_H__

namespace qucs {

class circuit;
class node;
class net;

namespace device {

  // creates external resistor circuit
  circuit *
    splitResistor (
      circuit * base, // calling circuit (this)
      circuit * res,  // additional resistor circuit (can be NULL)
      const char * c, // name of the additional circuit
      const char * n, // name of the inserted (internal) node
      int internal);  // number of new node (the original external node)

  // removes external resistor circuit
  void
    disableResistor (
      circuit * base, // calling circuit (this)
      circuit * res,  // additional resistor circuit
      int internal);  // number of new node (the original external node)

  // creates external capacitor circuit
  circuit *
    splitCapacitor (
      circuit * base, // calling circuit (this)
      circuit * cap,  // additional capacitor circuit (can be NULL)
      const char * c, // name of the additional circuit
      node * n1,      // first node of new capacitor
      node * n2);     // second node of new capacitor

  // removes external capacitor circuit
  void
    disableCapacitor (
      circuit * base, // calling circuit (this)
      circuit * cap); // additional capacitor circuit

  // checks whether circuit is enabled
  int
    deviceEnabled (
      circuit * c); // circuit to be checked

  // computes current and its derivative for a MOS pn-junction
  void
    pnJunctionMOS (
      double Upn, // pn-voltage
      double Iss, // saturation current
      double Ute, // temperature voltage
      double& I,  // result current
      double& g); // result derivative

  // computes current and its derivative for a bipolar pn-junction
  void
    pnJunctionBIP (
      double Upn, // pn-voltage
      double Iss, // saturation current
      double Ute, // temperature voltage
      double& I,  // result current
      double& g); // result derivative

  // limits the forward pn-voltage
  double
     pnVoltage (
       double Ud,     // current pn-voltage
       double Uold,   // previous pn-voltage
       double Ut,     // temperature voltage
       double Ucrit); // critical voltage

  // computes the exponential pn-junction current
  double
    pnCurrent (
      double Upn,  // pn-voltage
      double Iss,  // saturation current
      double Ute); // temperature voltage

  // computes the exponential pn-junction current's derivative
  double
    pnConductance (
      double Upn,  // pn-voltage
      double Iss,  // saturation current
      double Ute); // temperature voltage

  // computes pn-junction depletion capacitance
  double
    pnCapacitance (
      double Uj,  // pn-voltage
      double Cj,  // zero-bias capacitance
      double Vj,  // built-in potential
      double Mj,  // grading coefficient
      double Fc); // forward-bias coefficient

  // computes pn-junction depletion charge
  double
    pnCharge (
      double Uj,  // pn-voltage
      double Cj,  // zero-bias capacitance
      double Vj,  // built-in potential
      double Mj,  // grading coefficient
      double Fc); // forward-bias coefficient

  // computes pn-junction depletion capacitance
  double
    pnCapacitance (
      double Uj,  // pn-voltage
      double Cj,  // zero-bias capacitance
      double Vj,  // built-in potential
      double Mj); // grading coefficient

  // computes pn-junction depletion charge
  double
    pnCharge (
      double Uj,  // pn-voltage
      double Cj,  // zero-bias capacitance
      double Vj,  // built-in potential
      double Mj); // grading coefficient

  // compute critical voltage of pn-junction
  double
    pnCriticalVoltage (
      double Iss,  // saturation current
      double Ute); // temperature voltage

  // limits the forward fet-voltage
  double
    fetVoltage (
      double Ufet, // current fet voltage
      double Uold, // previous fet voltage
      double Uth); // threshold voltage

  // limits the drain-source voltage
  double
    fetVoltageDS (
      double Ufet,  // current fet voltage
      double Uold); // previous fet voltage

  // calculates the overlap capacitance for mosfet (meyer model)
  void
    fetCapacitanceMeyer (
      double Ugs,   // gate-source voltage
      double Ugd,   // gate-drain voltage
      double Uth,   // threshold voltage
      double Udsat, // drain-source saturation voltage
      double Phi,   // built-in potential
      double Cox,   // oxide capacitance
      double& Cgs,  // resulting gate-source capacitance
      double& Cgd,  // resulting gate-drain capacitance
      double& Cgb); // resulting gate-bulk capacitance

  // computes temperature dependency of energy bandgap
  double
    Egap (
      double T,            // temperature
      double Eg0 = Eg0Si); // bandgap at 0K

  // computes temperature dependency of intrinsic density
  double
    intrinsicDensity (
      double T,            // temperature
      double Eg0 = Eg0Si); // bandgap at 0K

  // calculates temperature dependence for saturation current
  double
    pnCurrent_T (
      double T1,       // reference temperature
      double T2,       // temperature
      double Is,       // saturation current
      double Eg,       // bandgap at 300K
      double N = 1,    // emission coefficient
      double Xti = 0); // temperature coefficient

  // calculates temperature dependence for junction potential
  double
    pnPotential_T (
      double T1,           // reference temperature
      double T2,           // temperature
      double Vj,           // built-in potential
      double Eg0 = Eg0Si); // bandgap at 0K

  // calculates temperature dependence for junction capacitance
  double
    pnCapacitance_T (
      double T1,  // reference temperature
      double T2,  // temperature
      double M,   // grading coefficient
      double VR,  // built-in potential ratio
      double Cj); // zero-bias capacitance

  // calculates temperature dependence for junction capacitance
  double
    pnCapacitance_F (
      double T1,  // reference temperature
      double T2,  // temperature
      double M,   // grading coefficient
      double VR); // built-in potential ratio: Vj(T2) / Vj(T1)

} // namespace device

} // namespace qucs

#endif /* __DEVICE_H__ */
