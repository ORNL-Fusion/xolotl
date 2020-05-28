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
constexpr double kBoltzmann = 8.61733240000000000E-5;

//! Pi, taken from "100000 digits of Pi,"
//! at http://www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html
constexpr double pi = 3.1415926535897932;

//! Lattice Parameter. Equal to 3.17 Angstroms, taken from Becquart et. al.
//! Journal of Nuclear Materials 403 (2010) 75–88. Given in units here of nm.
constexpr double tungstenLatticeConstant = 0.31700000000000000;

//! Lattice Parameter. Given in units here of nm.
constexpr double alloyLatticeConstant = 0.36000000000000000;

//! Lattice Parameter for UO2
constexpr double uraniumDioxydeLatticeConstant = 0.57400000000000000;

//! Lattice Parameter for Iron
constexpr double ironLatticeConstant = 0.28700000000000000;

//! Core radius. Given in units here of nm.
constexpr double alloyCoreRadius = 0.36000000000000000;

//! Single helium radius. Given in units here of nm.
constexpr double heliumRadius = 0.30000000000000000;

//! Single xenon radius. Given in units here of nm.
constexpr double xenonRadius = 0.30000000000000000;

// Tungsten heat coefficient = lambda / (rho * C) in nm2 s-1
constexpr double tungstenHeatCoefficient = 6.835e13;

// UO2 heat coefficient = lambda / (rho * C) in nm2 s-1
constexpr double uo2HeatCoefficient = 0.0;

// Iron heat coefficient = lambda / (rho * C) in nm2 s-1
constexpr double feHeatCoefficient = 0.0;

// Alloy heat coefficient = lambda / (rho * C) in nm2 s-1
constexpr double alloyHeatCoefficient = 0.0;

// Tungsten heat conductivity = lambda in W K-1 nm-1
constexpr double tungstenHeatConductivity = 173 * 1.0e-9;

// UO2 heat conductivity = lambda in W K-1 m-1
constexpr double uo2HeatConductivity = 0.0;

// Iron heat conductivity = lambda in W K-1 m-1
constexpr double feHeatConductivity = 0.0;

// Alloy heat conductivity = lambda in W K-1 m-1
constexpr double alloyHeatConductivity = 0.0;

// Burgers vector magnitudes for loops in Alloy case
// In lattice parameter units
constexpr double perfectBurgers = 0.5;
constexpr double faultedBurgers = 0.333;
constexpr double frankBurgers = 0.333;

// Sink strength for Alloy case in nm^-2
constexpr double alloysinkStrength = 1.0e-5;

} /* end namespace xolotlCore */
#endif /* CONSTANTS_H_ */
