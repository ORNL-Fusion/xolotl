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
    computeRate(TReactionNetwork& network)
    {
        return 1.0;
    }
};
}
}
