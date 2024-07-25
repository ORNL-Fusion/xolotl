#pragma once

#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace t91
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

	// Special case for V interacting with HeV
	if (lo1.isOnAxis(Species::V) or lo2.isOnAxis(Species::V)) {
		if (lo1[(int)Species::He] > 0 and lo1[(int)Species::V] > 0)
			rCore = ::xolotl::core::t91VCavRadius;
		else if (lo2[(int)Species::He] > 0 and lo2[(int)Species::V] > 0)
			rCore = ::xolotl::core::t91VCavRadius;
	}

	// Special case for I interacting with HeV
	if (lo1.isOnAxis(Species::I) or lo2.isOnAxis(Species::I)) {
		if (lo1[(int)Species::He] > 0 and lo1[(int)Species::V] > 0)
			rCore = ::xolotl::core::t91ICavRadius;
		else if (lo2[(int)Species::He] > 0 and lo2[(int)Species::V] > 0)
			rCore = ::xolotl::core::t91ICavRadius;
	}

	double kPlus = 4.0 * pi * (r0 + r1 + rCore) * (dc0 + dc1);

	return kPlus;
}
} // namespace t91

KOKKOS_INLINE_FUNCTION
double
T91ProductionReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_reactants[0]);
	auto cl1 = this->_clusterData->getCluster(_reactants[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return t91::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1,
		this->_clusterData->latticeParameter());
}

KOKKOS_INLINE_FUNCTION
void
T91ProductionReaction::computeFlux(
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
		double thermalVConc = exp(::xolotl::core::t91FormationEntropy) *
			exp(-::xolotl::core::t91FormationEnergy /
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
T91ProductionReaction::computePartialDerivatives(
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
		double thermalVConc = exp(::xolotl::core::t91FormationEntropy) *
			exp(-::xolotl::core::t91FormationEnergy /
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
T91ProductionReaction::computeReducedPartialDerivatives(
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
		double thermalVConc = exp(::xolotl::core::t91FormationEntropy) *
			exp(-::xolotl::core::t91FormationEnergy /
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
T91DissociationReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_products[0]);
	auto cl1 = this->_clusterData->getCluster(_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return t91::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1,
		this->_clusterData->latticeParameter());
}

KOKKOS_INLINE_FUNCTION
double
T91DissociationReaction::computeBindingEnergy(double time)
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	constexpr double heBinding[5] = {0.0, 0.0, 0.43, 0.95, 0.98};
	constexpr double vBinding[5] = {0.0, 0.0, 0.30, 0.37, 0.62};

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
		// He
		if (comp.isOnAxis(Species::He)) {
			if (prod1Comp.isOnAxis(Species::He) ||
				prod2Comp.isOnAxis(Species::He)) {
				be = heBinding[comp[Species::He]];
			}
		}
		// V
		else if (comp.isOnAxis(Species::V)) {
			auto size = comp[Species::V];
			if (size < 5)
				be = vBinding[size];
			else
				be = 2.02 -
					2.93 *
						(pow((double)size, 2.0 / 3.0) -
							pow((double)size - 1.0, 2.0 / 3.0));
		}
		// No I for now because only single I should be included in the
		// simulation HeV
		else if (comp[Species::He] > 0 and comp[Species::V] > 0) {
			auto amtHe = comp[Species::He], amtV = comp[Species::V];
			auto ratio = (double)amtHe / (double)amtV;
			if (prod1Comp.isOnAxis(Species::V) ||
				prod2Comp.isOnAxis(Species::V)) {
				be = 2.02 -
					2.93 *
						(pow((double)amtV, 2.0 / 3.0) -
							pow((double)amtV - 1.0, 2.0 / 3.0)) -
					0.1 * ratio * ratio + 1.59 * ratio;
			}
			if (prod1Comp.isOnAxis(Species::I) ||
				prod2Comp.isOnAxis(Species::I)) {
				be = 3.77 +
					2.93 *
						(pow((double)amtV, 2.0 / 3.0) -
							pow((double)amtV - 1.0, 2.0 / 3.0)) +
					0.09 * ratio * ratio - 1.35 * ratio;
			}
			if (prod1Comp.isOnAxis(Species::He) ||
				prod2Comp.isOnAxis(Species::He)) {
				be = 4.44 -
					1.99 *
						(pow((double)amtV, 2.0 / 3.0) -
							pow((double)amtV - 1.0, 2.0 / 3.0)) +
					0.07 * ratio * ratio - 1.06 * ratio;
			}
		}
	}
	// Grouped, only HeV should be grouped
	else {
		Composition lo = clReg.getOrigin();
		Composition hi = clReg.getUpperLimitPoint();
		Composition prod1Comp = prod1Reg.getOrigin();
		Composition prod2Comp = prod2Reg.getOrigin();
		// HeV
		double amtHe = (double)(lo[Species::He] + hi[Species::He] - 1) / 2.0,
			   amtV = (double)(lo[Species::V] + hi[Species::V] - 1) / 2.0;
		auto ratio = (double)amtHe / (double)amtV;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 2.02 -
				2.93 *
					(pow((double)amtV, 2.0 / 3.0) -
						pow((double)amtV - 1.0, 2.0 / 3.0)) -
				0.1 * ratio * ratio + 1.59 * ratio;
		}
		if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
			be = 3.77 +
				2.93 *
					(pow((double)amtV, 2.0 / 3.0) -
						pow((double)amtV - 1.0, 2.0 / 3.0)) +
				0.09 * ratio * ratio - 1.35 * ratio;
		}
		if (prod1Comp.isOnAxis(Species::He) ||
			prod2Comp.isOnAxis(Species::He)) {
			be = 4.44 -
				1.99 *
					(pow((double)amtV, 2.0 / 3.0) -
						pow((double)amtV - 1.0, 2.0 / 3.0)) +
				0.07 * ratio * ratio - 1.06 * ratio;
		}
	}

	return util::min(5.0, util::max(be, -5.0));
}

KOKKOS_INLINE_FUNCTION
double
T91SinkReaction::getSinkBias()
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
T91SinkReaction::getSinkStrength()
{
	auto bias = this->getSinkBias();

	double strength =
		::xolotl::core::t91DisloStrength * bias + ::xolotl::core::t91GBStrength;

	return strength;
}

KOKKOS_INLINE_FUNCTION
double
T91SinkReaction::computeRate(IndexType gridIndex, double time)
{
	auto cl = this->_clusterData->getCluster(_reactant);
	double dc = cl.getDiffusionCoefficient(gridIndex);

	double rate = this->getSinkStrength() * dc;

	return rate;
}
} // namespace network
} // namespace core
} // namespace xolotl
