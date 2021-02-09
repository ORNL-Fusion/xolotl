#pragma once

#include <xolotl/core/network/impl/NucleationReaction.tpp>
#include <xolotl/core/network/impl/ReSolutionReaction.tpp>
#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/core/network/impl/Reaction.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
KOKKOS_INLINE_FUNCTION
double
NEDissociationReaction::computeBindingEnergy()
{
	auto cl = this->_clusterData.getCluster(this->_reactant);
	auto prod1 = this->_clusterData.getCluster(this->_products[0]);
	auto prod2 = this->_clusterData.getCluster(this->_products[1]);
	return prod1.getFormationEnergy() + prod2.getFormationEnergy() -
		cl.getFormationEnergy();
}

KOKKOS_INLINE_FUNCTION
double
NESinkReaction::getSinkBias()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	double bias = 1.0;

	auto cl = this->_clusterData.getCluster(this->_reactant);

	auto clReg = cl.getRegion();
	if (clReg.isSimplex()) {
		Composition comp = clReg.getOrigin();
		if (comp.isOnAxis(Species::I)) {
			bias = 1.05;
		}
	}

	return bias;
}
} // namespace network
} // namespace core
} // namespace xolotl
