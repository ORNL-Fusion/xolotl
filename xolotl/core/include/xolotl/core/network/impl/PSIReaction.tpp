#pragma once

#include <xolotl/core/network/PSIClusterGenerator.h>
#include <xolotl/core/network/impl/Reaction.tpp>
#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/core/network/impl/TrapMutationReaction.tpp>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace psi
{
template <typename TRegion>
KOKKOS_INLINE_FUNCTION
double
getRate(const TRegion& pairCl0Reg, const TRegion& pairCl1Reg, const double r0,
	const double r1, const double dc0, const double dc1)
{
	constexpr double pi = ::xolotl::core::pi;

	double kPlus = 4.0 * pi * (r0 + r1) * (dc0 + dc1);

	return kPlus;
}
} // namespace psi

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSIProductionReaction<TSpeciesEnum>::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(this->_reactants[0]);
	auto cl1 = this->_clusterData->getCluster(this->_reactants[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return psi::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSIDissociationReaction<TSpeciesEnum>::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(this->_products[0]);
	auto cl1 = this->_clusterData->getCluster(this->_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return psi::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSIDissociationReaction<TSpeciesEnum>::computeBindingEnergy(double time)
{
	using psi::hasDeuterium;
	using psi::hasTritium;

	using NetworkType = typename Superclass::NetworkType;

	constexpr double beTableV1[10][7] = {
		// H:  1     2     3     4     5     6      // He:
		{0.0, 1.21, 1.17, 1.05, 0.93, 0.85, 0.60}, // 0
		{0.0, 1.00, 0.95, 0.90, 0.88, 0.80, 0.60}, // 1
		{0.0, 0.96, 0.92, 0.85, 0.84, 0.83, 0.50}, // 2
		{0.0, 0.86, 0.81, 0.69, 0.64, 0.65, 0.50}, // 3
		{0.0, 0.83, 0.80, 0.65, 0.60, 0.60, 0.55}, // 4
		{0.0, 0.83, 0.80, 0.60, 0.50, 0.50, 0.50}, // 5
		{0.0, 0.80, 0.70, 0.60, 0.50, 0.50, 0.50}, // 6
		{0.0, 0.80, 0.75, 0.65, 0.55, 0.55, 0.45}, // 7
		{0.0, 0.80, 0.80, 0.70, 0.65, 0.60, 0.55}, // 8
		{0.0, 0.80, 0.80, 0.75, 0.70, 0.65, 0.60}, // 9
	};

	constexpr double beTableV2[15][12] = {
		// H:  1     2     3     4     5     6     7     8     9     10    11 //
		// He:
		{0.0, 1.63, 1.31, 1.25, 1.16, 1.00, 1.00, 0.95, 0.95, 0.75, 0.70,
			0.65}, // 0
		{0.0, 1.30, 1.30, 1.24, 1.08, 0.95, 0.95, 0.95, 0.95, 0.75, 0.70,
			0.65}, // 1
		{0.0, 1.15, 1.14, 1.11, 1.14, 0.95, 0.95, 0.95, 0.90, 0.75, 0.70,
			0.65}, // 2
		{0.0, 1.12, 1.06, 0.99, 0.99, 0.90, 0.95, 0.90, 0.90, 0.70, 0.70,
			0.65}, // 3
		{0.0, 1.10, 1.06, 0.99, 0.99, 0.90, 0.95, 0.90, 0.90, 0.70, 0.65,
			0.65}, // 4
		{0.0, 1.10, 1.05, 0.99, 0.99, 0.90, 0.90, 0.90, 0.90, 0.70, 0.65,
			0.65}, // 5
		{0.0, 1.10, 1.05, 0.99, 0.99, 0.90, 0.90, 0.90, 0.85, 0.70, 0.65,
			0.60}, // 6
		{0.0, 1.05, 1.00, 0.95, 0.95, 0.90, 0.90, 0.90, 0.85, 0.65, 0.65,
			0.60}, // 7
		{0.0, 1.05, 1.00, 0.95, 0.95, 0.90, 0.90, 0.85, 0.85, 0.65, 0.65,
			0.60}, // 8
		{0.0, 1.05, 1.00, 0.95, 0.95, 0.85, 0.85, 0.85, 0.85, 0.65, 0.65,
			0.60}, // 9
		{0.0, 1.00, 0.95, 0.90, 0.90, 0.85, 0.85, 0.85, 0.80, 0.65, 0.60,
			0.60}, // 10
		{0.0, 0.95, 0.95, 0.90, 0.90, 0.85, 0.85, 0.85, 0.80, 0.65, 0.60,
			0.60}, // 11
		{0.0, 0.95, 0.90, 0.90, 0.85, 0.85, 0.85, 0.80, 0.80, 0.60, 0.60,
			0.55}, // 12
		{0.0, 0.90, 0.90, 0.85, 0.85, 0.85, 0.85, 0.80, 0.80, 0.60, 0.60,
			0.55}, // 13
		{0.0, 0.90, 0.90, 0.85, 0.85, 0.80, 0.80, 0.80, 0.70, 0.60, 0.60,
			0.55}, // 14
	};

	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;

	double be = 0.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);
	auto prod1 = this->_clusterData->getCluster(this->_products[0]);
	auto prod2 = this->_clusterData->getCluster(this->_products[1]);

	auto clReg = cl.getRegion();
	auto prod1Reg = prod1.getRegion();
	auto prod2Reg = prod2.getRegion();
	bool useTable = false;
	if (clReg.isSimplex()) {
		if (prod1Reg.isSimplex()) {
			auto orig1 = prod1Reg.getOrigin();
			if constexpr (hasDeuterium<Species> && hasTritium<Species>) {
				if (orig1.isOnAxis(Species::D) || orig1.isOnAxis(Species::T)) {
					useTable = true;
				}
			}
		}
		if (prod2Reg.isSimplex()) {
			auto orig2 = prod2Reg.getOrigin();
			if constexpr (hasDeuterium<Species> && hasTritium<Species>) {
				if (orig2.isOnAxis(Species::D) || orig2.isOnAxis(Species::T)) {
					useTable = true;
				}
			}
		}
	}

	if constexpr (hasDeuterium<Species> && hasTritium<Species>) {
		if (useTable) {
			Composition comp(clReg.getOrigin());
			auto hAmount = comp[Species::D] + comp[Species::T];
			if (comp[Species::V] == 1) {
				be = beTableV1[comp[Species::He]][hAmount];
			}
			else if (comp[Species::V] == 2) {
				be = beTableV2[comp[Species::He]][hAmount];
			}
		}
	}

	if (be == 0.0) {
		// Special case for V
		auto orig1 = prod1Reg.getOrigin();
		auto orig2 = prod2Reg.getOrigin();
		Composition comp(clReg.getOrigin());
		AmountType lowerV = 16, higherV = 31;
		AmountType minV = 1;
		for (auto i = 1; i < higherV; i++) {
			auto maxHe = psi::getMaxHePerV(i, 4.0);
			if (comp[Species::He] > maxHe)
				minV = i;
		}
		lowerV = util::max(lowerV, minV + 2);
		if ((orig1.isOnAxis(Species::V) || orig2.isOnAxis(Species::V)) &&
			(comp[Species::V] >= lowerV && comp[Species::V] <= higherV)) {
			// Get the be at 16 and 30
			Composition HeVComp(clReg.getOrigin());
			HeVComp[Species::V] = lowerV;
			auto fe1 = PSIClusterGenerator<TSpeciesEnum>::getHeVFormationEnergy(
				HeVComp);
			HeVComp[Species::V] = lowerV - 1;
			auto fe2 = PSIClusterGenerator<TSpeciesEnum>::getHeVFormationEnergy(
				HeVComp);
			Composition vComp{};
			vComp[Species::V] = 1;
			auto fe3 =
				PSIClusterGenerator<TSpeciesEnum>::getHeVFormationEnergy(vComp);
			auto be1 = fe2 + fe3 - fe1;
			HeVComp[Species::V] = higherV;
			fe1 = PSIClusterGenerator<TSpeciesEnum>::getHeVFormationEnergy(
				HeVComp);
			HeVComp[Species::V] = higherV - 1;
			fe2 = PSIClusterGenerator<TSpeciesEnum>::getHeVFormationEnergy(
				HeVComp);
			auto be2 = fe2 + fe3 - fe1;
			if (higherV - lowerV < 4)
				be = be2;
			else
				be = be1 +
					(comp[Species::V] - lowerV) * (be2 - be1) /
						(higherV - lowerV);
		}
		else {
			be = prod1.getFormationEnergy() + prod2.getFormationEnergy() -
				cl.getFormationEnergy();
		}
	}

	return util::max(be, -5.0);
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSISinkReaction<TSpeciesEnum>::getSinkBias()
{
	return 1.0;
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSISinkReaction<TSpeciesEnum>::getSinkStrength()
{
	constexpr double pi = ::xolotl::core::pi;
	double grainSize = 50000.0; // 50 um

	return 1.0 / (pi * grainSize * grainSize);
}
} // namespace network
} // namespace core
} // namespace xolotl
