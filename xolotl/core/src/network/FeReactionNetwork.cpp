#include <xolotl/core/network/FeReactionNetwork.h>
#include <xolotl/core/network/FeReactionNetwork.inl>

namespace xolotl
{
namespace core
{
namespace network
{
template
ReactionNetwork<FeReactionNetwork>::ReactionNetwork(
    const std::vector<AmountType>& maxSpeciesAmounts,
    const std::vector<SubdivisionRatio>& subdivisionRatios, IndexType gridSize,
    const options::IOptions& opts);

template
ReactionNetwork<FeReactionNetwork>::ReactionNetwork(
    const std::vector<AmountType>& maxSpeciesAmounts, IndexType gridSize,
    const options::IOptions& opts);

template
double
ReactionNetwork<FeReactionNetwork>::getTotalConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize);

template
double
ReactionNetwork<FeReactionNetwork>::getTotalRadiusConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize);

template
double
ReactionNetwork<FeReactionNetwork>::getTotalAtomConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize);

template
double
ReactionNetwork<FeReactionNetwork>::getTotalTrappedAtomConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize);
}
}
}
