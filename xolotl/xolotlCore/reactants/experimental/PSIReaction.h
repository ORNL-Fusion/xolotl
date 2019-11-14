#pragma once

#include <experimental/Reaction.h>

namespace xolotlCore
{
namespace experimental
{
class PSIReaction : public Reaction<PSIReaction>
{
public:
    using Reaction<PSIReaction>::Reaction;

    template <typename TReactionNetwork>
    double
    computeBindingEnergy(TReactionNetwork& network)
    {
        // static constexpr double bindingEnergyMap[2];

        return 1.0;
    }
};
}
}
