#pragma once

#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/core/network/impl/TransformReaction.tpp>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TRegion>
KOKKOS_INLINE_FUNCTION
double
getRate(const TRegion& pairCl0Reg, const TRegion& pairCl1Reg, const double r0,
	const double r1, const double dc0, const double dc1)
{
	constexpr double pi = ::xolotl::core::pi;
	const double zs = 4.0 * pi * (r0 + r1 + ::xolotl::core::fecrCoreRadius);

	using Species = typename TRegion::EnumIndex;
	auto lo0 = pairCl0Reg.getOrigin();
	auto lo1 = pairCl1Reg.getOrigin();

	if (lo0.isOnAxis(Species::I) and lo1.isOnAxis(Species::I))
		return 1.1 * zs * (dc0 + dc1);

	// react_extension
	if ((lo0.isOnAxis(Species::I) and lo1.isOnAxis(Species::V)) or
		(lo0.isOnAxis(Species::V) and lo1.isOnAxis(Species::I)))
		return 4.0 * pi * (r0 + r1 + ::xolotl::core::fecrCoreRadius + 0.25) *
			(dc0 + dc1);

	bool cl0IsSphere = (lo0.isOnAxis(Species::V) ||
			 lo0.isOnAxis(Species::Complex) || lo0.isOnAxis(Species::I)),
		 cl1IsSphere = (lo1.isOnAxis(Species::V) ||
			 lo1.isOnAxis(Species::Complex) || lo1.isOnAxis(Species::I));
	bool cl0IsTrap = lo0.isOnAxis(Species::Trap),
		 cl1IsTrap = lo1.isOnAxis(Species::Trap);

	// Simple case
	if (cl0IsSphere && cl1IsSphere) {
		return zs * (dc0 + dc1);
	}

	// With a trap
	if (cl0IsTrap || cl1IsTrap) {
		double rT = cl0IsTrap ? r0 : r1;
		double rL = cl0IsTrap ? r1 : r0;
		double sigma = pi * (r0 + r1) * (r0 + r1);
		if (rT > rL)
			sigma = 4.0 * pi * r0 * r1;
		return 2.0 * (dc0 + dc1) * sigma;
	}

	double rS = cl0IsSphere ? r0 : r1;
	double rL = cl0IsSphere ? r1 : r0;
	double rC = ::xolotl::core::fecrCoreRadius;
	if (!cl0IsSphere && !cl1IsSphere) {
		rS = 0.0;
		rL = r0 + r1;
		rC = ::xolotl::core::fecrCoalesceRadius;
	}
	double ratio = rL / (3.0 * (rS + rC));
	double p = 1.0 / (1.0 + ratio * ratio);
	double zl = 4.0 * pi * pi * rL / log(1.0 + 8.0 * rL / (rS + rC));
	double zd = -8.0 * pi * pi * rL /
		log(pi * ::xolotl::core::fecrSinkStrength * (rS + rC) * (rS + rC));

	double k_plus = (dc0 + dc1) * (p * zs + (1.0 - p) * std::max(zd, zl));
	double bias = 1.0;
	if (lo0.isOnAxis(Species::I) || lo1.isOnAxis(Species::I)) {
		if (lo0.isOnAxis(Species::Free) || lo1.isOnAxis(Species::Free))
			bias = 1.0;
		else if (lo0[static_cast<int>(Species::Loop)] > 0 ||
			lo1[static_cast<int>(Species::Loop)] > 0)
			bias = 1.05774;
		else
			bias = 1.05;
	}

	return k_plus * bias;
}

KOKKOS_INLINE_FUNCTION
double
FeCrProductionReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_reactants[0]);
	auto cl1 = this->_clusterData->getCluster(_reactants[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

KOKKOS_INLINE_FUNCTION
double
FeCrDissociationReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_products[0]);
	auto cl1 = this->_clusterData->getCluster(_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

KOKKOS_INLINE_FUNCTION
double
FeCrDissociationReaction::computeRate(IndexType gridIndex, double)
{
	double T = this->_clusterData->temperature(gridIndex);
	constexpr double pi = ::xolotl::core::pi;
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	double kPlus = this->asDerived()->getRateForProduction(gridIndex);
	double E_b = this->asDerived()->computeBindingEnergy();

	constexpr double k_B = ::xolotl::core::kBoltzmann;

	auto cl0 = this->_clusterData->getCluster(_products[0]);
	auto cl1 = this->_clusterData->getCluster(_products[1]);

	auto lo0 = cl0.getRegion().getOrigin();
	auto lo1 = cl1.getRegion().getOrigin();
	bool cl0IsSphere = (lo0.isOnAxis(Species::V) ||
			 lo0.isOnAxis(Species::Complex) || lo0.isOnAxis(Species::I)),
		 cl1IsSphere = (lo1.isOnAxis(Species::V) ||
			 lo1.isOnAxis(Species::Complex) || lo1.isOnAxis(Species::I));
	bool cl0IsTrap = lo0.isOnAxis(Species::Trap),
		 cl1IsTrap = lo1.isOnAxis(Species::Trap);

	double kMinus = kPlus * std::exp(-E_b / (k_B * T));

	// Recoil detrap
	if ((lo0.isOnAxis(Species::Trap) && lo1.isOnAxis(Species::Free)) ||
		(lo0.isOnAxis(Species::Free) && lo1.isOnAxis(Species::Trap))) {
		auto cl = this->_clusterData->getCluster(_reactant);
		Composition lo = cl.getRegion().getOrigin();
		double kRecoil = exp(-(double)lo[Species::Trapped] / 10000.0);
		if (lo[Species::Junction] > 0)
			kRecoil = exp(-(double)lo[Species::Junction] / 400.0);
		//			kRecoil = 0.0;

		return ((kMinus / this->_clusterData->atomicVolume()) +
				   ::xolotl::core::detrapFrequency) *
			kRecoil;
	}

	// Standard case
	if (cl0IsSphere and cl1IsSphere)
		return (kMinus / this->_clusterData->atomicVolume());

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();
	auto rCoal = ::xolotl::core::fecrCoalesceRadius;

	double a = r0 + r1 + rCoal;
	double b = std::max(r0 + r1 - rCoal, 0.0);
	double sigma = ::xolotl::core::pi * (a * a - b * b);

	const double jumpDistance = this->_clusterData->latticeParameter() * sqrt(3.0) / 2.0;

	return kMinus / (sigma * jumpDistance);

	// Trap case
	if (cl0IsTrap or cl1IsTrap) {
		double rT = cl0IsTrap ? r0 : r1;
		double rL = cl0IsTrap ? r1 : r0;
		sigma = pi * (r0 + r1) * (r0 + r1);
		if (rT > rL)
			sigma = 4.0 * pi * r0 * r1;
	}
	// Loop and sphere case
	else if (cl0IsSphere or cl1IsSphere) {
		double rA = cl0IsSphere ? r0 : r1;
		double rP = cl0IsSphere ? r1 : r0;
		double rO = rP + ::xolotl::core::fecrCoreRadius;
		for (auto q : _align4by4) {
			double aOut = rA + rO;
			double bOut = rA + rO * q;
			double aIn = std::max(rA - rO, 0.0);
			double bIn = std::max(rA - rO * q, 0.0);

			sigma += _q4by4 * pi * (aOut * bOut - aIn * bIn);
		}
	}
	// Loop and loop
	else {
		double rO = ::xolotl::core::fecrCoalesceRadius;
		for (auto q : _align4by4) {
			double aOut = r0 + r1 + rO;
			double bOut = r0 + r1 * q + rO;
			double aIn = std::max(r0 + r1 - rO, 0.0);
			double bIn = std::max(r0 + r1 * q - rO, 0.0);

			sigma += _q4by4 * pi * (aOut * bOut - aIn * bIn);
		}
	}

	return kMinus / (sigma * jumpDistance);
}

KOKKOS_INLINE_FUNCTION
double
FeCrDissociationReaction::computeBindingEnergy(double time)
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	double be = 5.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);
	auto prod1 = this->_clusterData->getCluster(this->_products[0]);
	auto prod2 = this->_clusterData->getCluster(this->_products[1]);

	auto clReg = cl.getRegion();
	auto prod1Reg = prod1.getRegion();
	auto prod2Reg = prod2.getRegion();
	Composition lo = clReg.getOrigin();
	Composition hi = clReg.getUpperLimitPoint();
	Composition prod1Comp = prod1Reg.getOrigin();
	Composition prod2Comp = prod2Reg.getOrigin();
	if (lo.isOnAxis(Species::V)) {
		double n = (double)(lo[Species::V] + hi[Species::V] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::V) && prod2Comp.isOnAxis(Species::V)) {
			be = 1.73 - 2.59 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
		}
	}
	else if (lo[Species::Complex] > 0) {
		be = 0.01;
	}
	else if (lo.isOnAxis(Species::I)) {
		double n = (double)(lo[Species::I] + hi[Species::I] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::I) && prod2Comp.isOnAxis(Species::I)) {
			be = 4.33 - 5.76 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
		}
	}
	else if (lo.isOnAxis(Species::Free)) {
		double n = (double)(lo[Species::Free] + hi[Species::Free] - 1) / 2.0;
		be = 4.33 - 5.76 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
	}
	else if (lo[Species::Trapped] > 0) {
		if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
			double n =
				(double)(lo[Species::Trapped] + hi[Species::Trapped] - 1) / 2.0;
			be = 4.33 - 5.76 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
		}
		else if (prod1Comp.isOnAxis(Species::Trap) ||
			prod2Comp.isOnAxis(Species::Trap)) {
			be = 1.2;
		}
	}
	else if (lo[Species::Junction] > 0) {
		be = 2.5;
	}

	return util::min(5.0, util::max(be, -5.0));
}

KOKKOS_INLINE_FUNCTION
double
FeCrSinkReaction::getSinkBias()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	double bias = 1.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);

	auto clReg = cl.getRegion();
	if (clReg.isSimplex()) {
		Composition comp = clReg.getOrigin();
		if (comp.isOnAxis(Species::I)) {
			bias = 1.05;
		}
	}

	return bias;
}

KOKKOS_INLINE_FUNCTION
double
FeCrSinkReaction::getSinkStrength()
{
	auto density = this->_clusterData->sinkDensity(); // nm-2
	auto portion = this->_clusterData->sinkPortion();
	auto r = 1.0 / sqrt(::xolotl::core::pi * density);
	auto rCore = ::xolotl::core::fecrCoreRadius;
	auto temperature = this->_clusterData->temperature(0);
	constexpr double K = 170.0e9; // GPa
	constexpr double nu = 0.29;
	constexpr double b = 0.25; // nm
	double deltaV = 1.67 * this->_clusterData->atomicVolume();
	constexpr double a0 = 0.91, a1 = -2.16, a2 = -0.92; // Random dipole
	//	constexpr double a0 = 0.87, a1 = -5.12, a2 = -0.77; // Full network

	double L = (K * b * deltaV * (1.0 - 2.0 * nu)) /
		(2.0 * ::xolotl::core::pi * ::xolotl::core::kBoltzmann * temperature *
			(1.0 - nu));

	double delta = std::sqrt(rCore * rCore + (L * L) / 4.0);

	double Z =
		(2.0 * ::xolotl::core::pi * (a0 + a1 * (rCore / r)) *
			(1.0 +
				portion *
					((std::log(r / rCore) *
						 (a0 * r + a1 * delta + a2 * (delta - rCore))) /
							(std::log(r / delta) * (a0 * r + a1 * rCore)) -
						1.0))) /
		(std::log(r / rCore));

	return density * Z;
}

KOKKOS_INLINE_FUNCTION
double
FeCrTransformReaction::getSize()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	auto cl = this->_clusterData->getCluster(this->_reactant);

	auto clReg = cl.getRegion();
	if (clReg.isSimplex()) {
		Composition comp = clReg.getOrigin();
		return comp[Species::Junction];
	}

	return 0.0;
}

KOKKOS_INLINE_FUNCTION
double
FeCrTransformReaction::getExponent()
{
	// Same for Loop and Trapped

	return 2.0;
}

KOKKOS_INLINE_FUNCTION
double
FeCrTransformReaction::getBarrier()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	auto cl = this->_clusterData->getCluster(this->_product);

	auto clReg = cl.getRegion();
	if (clReg.isSimplex()) {
		Composition comp1 = this->_clusterData->getCluster(this->_reactant)
								.getRegion()
								.getOrigin();
		Composition comp = clReg.getOrigin();
		if (comp[Species::Loop] > 0)
			return 0.8; // Loop
		return 0.75; // Trapped
	}

	return 0.0;
}
} // namespace network
} // namespace core
} // namespace xolotl
