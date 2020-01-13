#pragma once

#include <experimental/NETraits.h>

namespace xolotlCore
{
namespace experimental
{
class NEReactionNetwork;

class NEReaction : public Reaction<NEReactionNetwork, NEReaction>
{
public:
    using NetworkType = NEReactionNetwork;
    using Superclass = Reaction<NetworkType, NEReaction>;

    using Superclass::Superclass;

    KOKKOS_INLINE_FUNCTION
    double
    computeBindingEnergy()
    {
        auto cl = this->_clusterData.getCluster(this->_reactants[0]);
        auto prod1 = this->_clusterData.getCluster(this->_products[0]);
        auto prod2 = this->_clusterData.getCluster(this->_products[1]);
        return prod1.getFormationEnergy() + prod2.getFormationEnergy() -
            cl.getFormationEnergy();
    }
};
}
}
