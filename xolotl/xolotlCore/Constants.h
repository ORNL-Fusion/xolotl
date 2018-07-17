/*
 * Constants.h
 *
 *  Created on: May 6, 2013
 *      Author: bkj
 */

#include <memory>
#include <cmath>

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

namespace xolotlCore {

//! Definitions of fundamental constants used in Xolotl

//! The Boltzmann constant in units of eV K^-1.
static const double kBoltzmann = 8.61733240000000000E-5;

//! Pi, taken from "100000 digits of Pi,"
//! at http://www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html
static const double pi = 3.1415926535897932;

//! Lattice Parameter. Equal to 3.17 Angstroms, taken from Becquart et. al.
//! Journal of Nuclear Materials 403 (2010) 75–88. Given in units here of nm.
static const double tungstenLatticeConstant = 0.31700000000000000;

//! Lattice Parameter for UO2
static const double uraniumDioxydeLatticeConstant = 0.57400000000000000;

//! Lattice Parameter for Iron
static const double ironLatticeConstant = 0.28700000000000000;

// Tungsten heat coefficient = lambda / (rho * C) in nm2 s-1
static const double tungstenHeatCoefficient = 6.835e13;

// UO2 heat coefficient = lambda / (rho * C) in nm2 s-1
static const double uo2HeatCoefficient = 0.0;

// Iron heat coefficient = lambda / (rho * C) in nm2 s-1
static const double feHeatCoefficient = 0.0;

// Tungsten heat conductivity = lambda in W K-1 nm-1
static const double tungstenHeatConductivity = 173 * 1.0e-9;

// UO2 heat conductivity = lambda in W K-1 m-1
static const double uo2HeatConductivity = 0.0;

// Iron heat conductivity = lambda in W K-1 m-1
static const double feHeatConductivity = 0.0;

//! Parameters for biased sink in the iron case
static const double reactionRadius = ironLatticeConstant
		* pow((3.0) / pi, (1.0 / 3.0)) * 0.5;
static const double r0 = ironLatticeConstant * 0.75 * sqrt(3.0);
static const double rho = 0.0003;
static const double sinkStrength = -4.0 * pi * rho
		/ log(pi * rho * pow(reactionRadius + r0, 2.0));
static const double sinkBias = 1.05;

} /* end namespace xolotlCore */
#endif /* CONSTANTS_H_ */
