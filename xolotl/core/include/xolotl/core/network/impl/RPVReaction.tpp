#pragma once

#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace rpv
{
template <typename TRegion>
KOKKOS_INLINE_FUNCTION
double
getRate(const TRegion& pairCl0Reg, const TRegion& pairCl1Reg, const double r0,
	const double r1, const double dc0, const double dc1,
	const double latticeConst)
{
	constexpr double pi = ::xolotl::core::pi;
	double rCore = 0.0;

	auto lo1 = pairCl0Reg.getOrigin();
	auto lo2 = pairCl1Reg.getOrigin();

	using Species = typename TRegion::EnumIndex;

	// Recombination
	if (lo1.isOnAxis(Species::I) and lo2.isOnAxis(Species::V)) {
		if (lo1[(int)Species::I] == 1 and lo2[(int)Species::V] == 1)
			return 500.0 * latticeConst * (dc0 + dc1);
	}
	if (lo2.isOnAxis(Species::I) and lo1.isOnAxis(Species::V)) {
		if (lo2[(int)Species::I] == 1 and lo1[(int)Species::V] == 1)
			return 500.0 * latticeConst * (dc0 + dc1);
	}

	double kPlus = 4.0 * pi * (r0 + r1 + rCore) * (dc0 + dc1);

	return kPlus;
}
} // namespace rpv

KOKKOS_INLINE_FUNCTION
double
RPVProductionReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_reactants[0]);
	auto cl1 = this->_clusterData->getCluster(_reactants[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	auto rate = rpv::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1,
		this->_clusterData->latticeParameter());

	auto prod = this->_clusterData->getCluster(_products[0]);
	auto loProd = prod.getRegion().getOrigin();

	// Split the rate for loops
	auto lo1 = cl0.getRegion().getOrigin();
	auto lo2 = cl1.getRegion().getOrigin();
	if (lo1[(int)Species::I] > 3 and lo2[(int)Species::I] > 3) {
		auto prod = this->_clusterData->getCluster(_products[0]);
		if (prod.getRegion().getOrigin().isOnAxis(Species::I)) {
			// Check the sizes
			auto loProd = prod.getRegion().getOrigin();
			if (loProd[(int)Species::I] >= 40) {
				return rate * 0.5;
			}
		}
		if (prod.getRegion().getOrigin().isOnAxis(Species::Loop)) {
			// Check the sizes
			auto loProd = prod.getRegion().getOrigin();
			if (loProd[(int)Species::Loop] >= 40) {
				return rate * 0.5;
			}
		}
	}

	return rate;
}

KOKKOS_INLINE_FUNCTION
void
RPVProductionReaction::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	int nProd = 0;
	for (auto prodId : _products) {
		if (prodId != invalidIndex) {
			++nProd;
		}
	}

	if (nProd == 0) {
		// Compute thermal vacancy concentration
		double omega = this->_clusterData->atomicVolume();
		double thermalVConc = exp(::xolotl::core::feBCCFormationEntropy) *
			exp(-::xolotl::core::feBCCFormationEnergy /
				(::xolotl::core::kBoltzmann *
					this->_clusterData->temperature(gridIndex))) /
			omega;

		// Which one is V?
		auto cl0 = this->_clusterData->getCluster(_reactants[0]);
		auto cl1 = this->_clusterData->getCluster(_reactants[1]);
		Composition cl0Comp = cl0.getRegion().getOrigin();
		Composition cl1Comp = cl1.getRegion().getOrigin();
		auto cV = cl0Comp.isOnAxis(Species::V) ? concentrations[_reactants[0]] :
												 concentrations[_reactants[1]];
		auto cI = cl0Comp.isOnAxis(Species::V) ? concentrations[_reactants[1]] :
												 concentrations[_reactants[0]];

		double f = this->_rate(gridIndex) * cI * (cV + thermalVConc);

		Kokkos::atomic_sub(&fluxes[_reactants[0]], f);
		Kokkos::atomic_sub(&fluxes[_reactants[1]], f);

		return;
	}

	ProductionReaction::computeFlux(concentrations, fluxes, gridIndex);
}

KOKKOS_INLINE_FUNCTION
void
RPVProductionReaction::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	int nProd = 0;
	for (auto prodId : _products) {
		if (prodId != invalidIndex) {
			++nProd;
		}
	}

	if (nProd == 0) {
		// Compute thermal vacancy concentration
		double omega = this->_clusterData->atomicVolume();
		double thermalVConc = exp(::xolotl::core::feBCCFormationEntropy) *
			exp(-::xolotl::core::feBCCFormationEnergy /
				(::xolotl::core::kBoltzmann *
					this->_clusterData->temperature(gridIndex))) /
			omega;

		// Which one is V?
		auto cl0 = this->_clusterData->getCluster(_reactants[0]);
		auto cl1 = this->_clusterData->getCluster(_reactants[1]);
		Composition cl0Comp = cl0.getRegion().getOrigin();
		Composition cl1Comp = cl1.getRegion().getOrigin();
		auto cV = cl0Comp.isOnAxis(Species::V) ? concentrations[_reactants[0]] :
												 concentrations[_reactants[1]];
		auto cI = cl0Comp.isOnAxis(Species::V) ? concentrations[_reactants[1]] :
												 concentrations[_reactants[0]];

		double f = this->_rate(gridIndex);

		if (cl0Comp.isOnAxis(Species::V)) {
			// First partial (V)
			Kokkos::atomic_sub(&values(_connEntries[0][0][0][0]), f * cI);
			Kokkos::atomic_sub(&values(_connEntries[1][0][0][0]), f * cI);

			// Second partial (I)
			Kokkos::atomic_sub(
				&values(_connEntries[0][0][1][0]), f * (cV + thermalVConc));
			Kokkos::atomic_sub(
				&values(_connEntries[1][0][1][0]), f * (cV + thermalVConc));
		}
		else {
			// First partial (I)
			Kokkos::atomic_sub(
				&values(_connEntries[0][0][0][0]), f * (cV + thermalVConc));
			Kokkos::atomic_sub(
				&values(_connEntries[1][0][0][0]), f * (cV + thermalVConc));

			// Second partial (V)
			Kokkos::atomic_sub(&values(_connEntries[0][0][1][0]), f * cI);
			Kokkos::atomic_sub(&values(_connEntries[1][0][1][0]), f * cI);
		}

		return;
	}

	ProductionReaction::computePartialDerivatives(
		concentrations, values, gridIndex);
}

KOKKOS_INLINE_FUNCTION
void
RPVProductionReaction::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	int nProd = 0;
	for (auto prodId : _products) {
		if (prodId != invalidIndex) {
			++nProd;
		}
	}

	if (nProd == 0) {
		// Compute thermal vacancy concentration
		double omega = this->_clusterData->atomicVolume();
		double thermalVConc = exp(::xolotl::core::feBCCFormationEntropy) *
			exp(-::xolotl::core::feBCCFormationEnergy /
				(::xolotl::core::kBoltzmann *
					this->_clusterData->temperature(gridIndex))) /
			omega;

		// Which one is V?
		auto cl0 = this->_clusterData->getCluster(_reactants[0]);
		auto cl1 = this->_clusterData->getCluster(_reactants[1]);
		Composition cl0Comp = cl0.getRegion().getOrigin();
		Composition cl1Comp = cl1.getRegion().getOrigin();
		auto cV = cl0Comp.isOnAxis(Species::V) ? concentrations[_reactants[0]] :
												 concentrations[_reactants[1]];
		auto cI = cl0Comp.isOnAxis(Species::V) ? concentrations[_reactants[1]] :
												 concentrations[_reactants[0]];

		double f = this->_rate(gridIndex);

		if (cl0Comp.isOnAxis(Species::V)) {
			// First partial (V)
			Kokkos::atomic_sub(&values(_connEntries[0][0][0][0]), f * cI);

			// Second partial (I)
			Kokkos::atomic_sub(
				&values(_connEntries[1][0][1][0]), f * (cV + thermalVConc));
		}
		else {
			// First partial (I)
			Kokkos::atomic_sub(
				&values(_connEntries[0][0][0][0]), f * (cV + thermalVConc));

			// Second partial (V)
			Kokkos::atomic_sub(&values(_connEntries[1][0][1][0]), f * cI);
		}

		return;
	}

	ProductionReaction::computeReducedPartialDerivatives(
		concentrations, values, gridIndex);
}

KOKKOS_INLINE_FUNCTION
double
RPVDissociationReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_products[0]);
	auto cl1 = this->_clusterData->getCluster(_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return rpv::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1,
		this->_clusterData->latticeParameter());
}

KOKKOS_INLINE_FUNCTION
double
RPVDissociationReaction::computeBindingEnergy(double time)
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	constexpr double vBinding[5] = {0.0, 0.0, 0.30, 0.37, 0.62};

	double be = 5.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);

	auto clReg = cl.getRegion();
	Composition lo = clReg.getOrigin();
	Composition hi = clReg.getUpperLimitPoint();

	// V
	if (lo.isOnAxis(Species::V)) {
		double n = (double)(lo[Species::V] + hi[Species::V] - 1) * 0.5;

		if (n < 5)
			be = vBinding[(IndexType)n];
		else
			be = 2.02 - 2.93 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
	}
	// No I for now because only single I

	return util::min(5.0, util::max(be, -5.0));
}

KOKKOS_INLINE_FUNCTION
double
RPVSinkReaction::getSinkBias()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	double bias = 1.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);

	auto clReg = cl.getRegion();
	if (clReg.isSimplex()) {
		Composition comp = clReg.getOrigin();
		if (comp.isOnAxis(Species::I)) {
			bias = 1.1;
		}
	}

	return bias;
}

KOKKOS_INLINE_FUNCTION
double
RPVSinkReaction::getSinkStrength()
{
	auto bias = this->getSinkBias();

	double strength = ::xolotl::core::feBCCDisloStrength * bias +
		::xolotl::core::feBCCGBStrength;

	return strength;
}

KOKKOS_INLINE_FUNCTION
double
RPVSinkReaction::computeRate(IndexType gridIndex, double time)
{
	auto cl = this->_clusterData->getCluster(_reactant);
	double dc = cl.getDiffusionCoefficient(gridIndex);

	double rate = this->getSinkStrength() * dc;

	return rate;
}
} // namespace network
} // namespace core
} // namespace xolotl
