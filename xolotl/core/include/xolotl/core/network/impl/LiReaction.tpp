#pragma once

#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace li
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
} // namespace li

KOKKOS_INLINE_FUNCTION
double
LiProductionReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_reactants[0]);
	auto cl1 = this->_clusterData->getCluster(_reactants[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return li::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

KOKKOS_INLINE_FUNCTION
double
LiDissociationReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_products[0]);
	auto cl1 = this->_clusterData->getCluster(_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return li::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

KOKKOS_INLINE_FUNCTION
double
LiDissociationReaction::computeBindingEnergy(double time)
{
        return 1.6; // Dissociation energy
	double be = 5.0;

	return util::min(5.0, util::max(be, -5.0));
}
} // namespace network
} // namespace core
} // namespace xolotl
