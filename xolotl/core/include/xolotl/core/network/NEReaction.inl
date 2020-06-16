#pragma once

#include <xolotl/core/network/ReSolutionReaction.inl>

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
}
}
}
