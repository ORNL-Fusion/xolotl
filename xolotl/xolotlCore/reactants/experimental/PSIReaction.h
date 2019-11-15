#pragma once

namespace xolotlCore
{
namespace experimental
{
template <typename TSpeciesEnum>
class PSIReaction : public
    ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>
        ::template Reaction<PSIReaction<TSpeciesEnum>>
{
public:
    using NetworkType = PSIReactionNetwork<TSpeciesEnum>;

    using Superclass = typename ReactionNetwork<NetworkType>
        ::template Reaction<PSIReaction<TSpeciesEnum>>;
    using Superclass::Superclass;

    double
    computeBindingEnergy()
    {
        // static constexpr double bindingEnergyMap[2];

        return 1.0;
    }
};
}
}
