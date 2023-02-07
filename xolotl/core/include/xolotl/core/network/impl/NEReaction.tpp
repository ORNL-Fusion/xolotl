#pragma once

#include <xolotl/core/network/impl/FullReSolutionReaction.tpp>
#include <xolotl/core/network/impl/NucleationReaction.tpp>
#include <xolotl/core/network/impl/PartialReSolutionReaction.tpp>
#include <xolotl/core/network/impl/Reaction.tpp>

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

	double kPlus = 4.0 * pi * (r0 + r1) * (dc0 + dc1);

	return kPlus;
}

KOKKOS_INLINE_FUNCTION
double
NEProductionReaction::getRateForProduction(IndexType gridIndex)
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
NEDissociationReaction::getRateForProduction(IndexType gridIndex)
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
NEDissociationReaction::computeBindingEnergy(double time)
{
	auto cl = this->_clusterData->getCluster(this->_reactant);
	auto prod1 = this->_clusterData->getCluster(this->_products[0]);
	auto prod2 = this->_clusterData->getCluster(this->_products[1]);
	return prod1.getFormationEnergy() + prod2.getFormationEnergy() -
		cl.getFormationEnergy();
}

KOKKOS_INLINE_FUNCTION
void
NEFullReSolutionReaction::setSize()
{
	auto cl = this->_clusterData->getCluster(this->_reactant);
	auto reg = cl.getRegion();
	Composition lo = reg.getOrigin();
	Composition hi = reg.getUpperLimitPoint();

	this->_size =
		((double)(hi[Species::Xe] - 1) * (double)hi[Species::Xe] / 2.0 -
			(double)(lo[Species::Xe] - 1) * (double)lo[Species::Xe] / 2.0) /
		reg.volume();
}
} // namespace network
} // namespace core
} // namespace xolotl
