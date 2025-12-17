#pragma once

#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/core/network/impl/TrapReaction.tpp>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace fe
{
template <typename TRegion>
KOKKOS_INLINE_FUNCTION
double
getRate(const TRegion& pairCl0Reg, const TRegion& pairCl1Reg, const double r0,
	const double r1, const double dc0, const double dc1)
{
	constexpr double pi = ::xolotl::core::pi;

	double kPlus =
		4.0 * pi * (r0 + r1 + ::xolotl::core::feCrCoreRadius) * (dc0 + dc1);

	return kPlus;
}
} // namespace fe

KOKKOS_INLINE_FUNCTION
double
FeProductionReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_reactants[0]);
	auto cl1 = this->_clusterData->getCluster(_reactants[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return fe::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

KOKKOS_INLINE_FUNCTION
double
FeDissociationReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_products[0]);
	auto cl1 = this->_clusterData->getCluster(_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return fe::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

KOKKOS_INLINE_FUNCTION
double
FeDissociationReaction::computeBindingEnergy(double time)
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	constexpr double heTrapTable[9] = {
		0.0, 4.31, 2.90, 2.02, 1.09, 0.58, 0.13, -0.25, -0.59};

	constexpr double heHeTable[9] = {
		0.0, 0.0, 0.43, 0.95, 0.98, 1.0, 1.0, 1.0, 1.0};

	constexpr double vVTable[5] = {0.0, 0.0, 0.30, 0.37, 0.62};

	constexpr double bubbleHeTable[4][4] = {{2.3, 2.84, 3.3, 3.84},
		{1.84, 2.75, 2.96, 3.12}, {1.83, 2.07, 2.91, 3.16},
		{1.91, 2.36, 2.57, 3.05}};

	constexpr double bubbleV1Table[8] = {
		3.71, 4.59, 5.52, 6.02, 6.48, 6.86, 7.20, 7.49};

	constexpr double bubbleVTable[3][4] = {{0.78, 1.61, 1.85, 2.30},
		{0.83, 1.04, 1.80, 2.03}, {1.16, 1.32, 1.57, 1.97}};

	constexpr double bubbleITable[3][4] = {{5.83, 5.0, 4.76, 4.31},
		{5.78, 5.57, 4.81, 4.58}, {5.45, 5.29, 5.04, 4.64}};

	double be = 5.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);
	auto prod1 = this->_clusterData->getCluster(this->_products[0]);
	auto prod2 = this->_clusterData->getCluster(this->_products[1]);

	auto clReg = cl.getRegion();
	auto prod1Reg = prod1.getRegion();
	auto prod2Reg = prod2.getRegion();
	if (clReg.isSimplex() && prod1Reg.isSimplex() && prod2Reg.isSimplex()) {
		Composition comp = clReg.getOrigin();
		Composition prod1Comp = prod1Reg.getOrigin();
		Composition prod2Comp = prod2Reg.getOrigin();
		if (comp.isOnAxis(Species::He)) {
			if (prod1Comp.isOnAxis(Species::He) ||
				prod2Comp.isOnAxis(Species::He)) {
				be = heHeTable[comp[Species::He]];
			}
			if (prod1Comp.isOnAxis(Species::I) ||
				prod2Comp.isOnAxis(Species::I)) {
				be = heTrapTable[comp[Species::He]];
			}
		}
		else if (comp.isOnAxis(Species::V)) {
			auto size = comp[Species::V];
			if (size < 5)
				be = vVTable[size];
			else {
				be = 1.73 -
					2.59 *
						(pow((double)size, 2.0 / 3.0) -
							pow((double)size - 1.0, 2.0 / 3.0));
			}
		}
		else if (comp.isOnAxis(Species::I)) {
			// Nothing
		}
		else {
			// HeV
			auto amtHe = comp[Species::He], amtV = comp[Species::V];
			double omega = this->_clusterData->atomicVolume();
			double T = this->_clusterData->temperature(0);
			constexpr double k_B = ::xolotl::core::kBoltzmann;
			double dg = 0.3135 * (0.8542 - 0.03996 * log(9.16 / T));
			double eta = ::xolotl::core::pi * dg * dg * dg * (double)amtHe /
				(6.0 * omega * (double)amtV);
			double z = (1.0 + eta + eta * eta * (1.0 - eta)) /
				((1.0 - eta) * (1.0 - eta) * (1.0 - eta));
			double p = (double)amtHe * z * k_B * T / (double)amtV;
			// HeV -> V
			if (prod1Comp.isOnAxis(Species::V) ||
				prod2Comp.isOnAxis(Species::V)) {
				return 5.0;
				if (amtV == 1 and amtHe < 9)
					be = bubbleV1Table[amtHe - 1];
				else if (amtV < 5 and amtHe < 5) {
					be = bubbleVTable[amtV - 2][amtHe - 1];
				}
				else
					be = 1.73 -
						2.59 *
							(pow((double)amtV, 2.0 / 3.0) -
								pow((double)amtV - 1.0, 2.0 / 3.0)) +
						p * omega;
			}
			// HeV -> I
			if (prod1Comp.isOnAxis(Species::I) ||
				prod2Comp.isOnAxis(Species::I)) {
				if (amtV < 4 and amtHe < 5)
					be = bubbleITable[amtV - 1][amtHe - 1];
				else
					be = 4.88 +
						2.59 *
							(pow((double)amtV, 2.0 / 3.0) -
								pow((double)amtV - 1.0, 2.0 / 3.0)) -
						p * omega;
			}
		}
	}
	else {
		Composition lo = clReg.getOrigin();
		Composition hi = clReg.getUpperLimitPoint();
		Composition prod1Comp = prod1Reg.getOrigin();
		Composition prod2Comp = prod2Reg.getOrigin();
		// HeV
		auto amtHe = (double)(lo[Species::He] + hi[Species::He] - 1) / 2.0;
		auto amtV = (double)(lo[Species::V] + hi[Species::V] - 1) / 2.0;
		double omega = this->_clusterData->atomicVolume();
		double T = this->_clusterData->temperature(0);
		constexpr double k_B = ::xolotl::core::kBoltzmann;
		double dg = 0.3135 * (0.8542 - 0.03996 * log(9.16 / T));
		double eta = ::xolotl::core::pi * dg * dg * dg * (double)amtHe /
			(6.0 * omega * (double)amtV);
		double z = (1.0 + eta + eta * eta * (1.0 - eta)) /
			((1.0 - eta) * (1.0 - eta) * (1.0 - eta));
		double p = (double)amtHe * z * k_B * T / (double)amtV;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			return 5.0;
			be = 1.73 -
				2.59 * (pow(amtV, 2.0 / 3.0) - pow(amtV - 1.0, 2.0 / 3.0)) +
				p * omega;
		}
		if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
			be = 4.88 +
				2.59 * (pow(amtV, 2.0 / 3.0) - pow(amtV - 1.0, 2.0 / 3.0)) -
				p * omega;
		}
	}

	return util::min(5.0, util::max(be, -5.0));
}

KOKKOS_INLINE_FUNCTION
double
FeSinkReaction::getSinkBias()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	double bias = 1.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);

	auto clReg = cl.getRegion();
	if (clReg.isSimplex()) {
		Composition comp = clReg.getOrigin();
		if (comp.isOnAxis(Species::I)) {
			bias = 1.2;
		}
		if (comp.isOnAxis(Species::He)) {
			bias = 0.8;
		}
	}

	return bias;
}

KOKKOS_INLINE_FUNCTION
double
FeSinkReaction::getSinkStrength()
{
	auto cl = this->_clusterData->getCluster(this->_reactant);
	double r = cl.getReactionRadius();
	double latticeParameter = this->_clusterData->latticeParameter();
	double r0 = latticeParameter * 0.75 * sqrt(3.0);
	double rho = 0.00025;
	constexpr double pi = ::xolotl::core::pi;

	double strength =
		-4.0 * pi * rho * (r + r0) / log(pi * rho * (r + r0) * (r + r0));

	return strength;
}

KOKKOS_INLINE_FUNCTION
IdType
FeTrapReaction::getId()
{
	using Composition = typename Superclass::Composition;

	// Here the id is the size of the trap
	auto cl = this->_clusterData->getCluster(this->_trapped);
	Composition comp = cl.getRegion().getOrigin();
	return comp[Species::Trap] - 1;
}

KOKKOS_INLINE_FUNCTION
double
FeTrapReaction::getStrength()
{
	return strength;
}

KOKKOS_INLINE_FUNCTION
double
FeTrapReaction::getEnergy()
{
	return energy;
}

KOKKOS_INLINE_FUNCTION
double
FeTrapReaction::getFrequency()
{
	return frequency;
}

KOKKOS_INLINE_FUNCTION
double
FeTrapReaction::getDensity()
{
	return density;
}
} // namespace network
} // namespace core
} // namespace xolotl
