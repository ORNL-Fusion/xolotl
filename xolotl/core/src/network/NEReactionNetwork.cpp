#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/core/network/impl/NEReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template
ReactionNetwork<NEReactionNetwork>::ReactionNetwork(
    const std::vector<AmountType>& maxSpeciesAmounts,
    const std::vector<SubdivisionRatio>& subdivisionRatios, IndexType gridSize,
    const options::IOptions& opts);

template
ReactionNetwork<NEReactionNetwork>::ReactionNetwork(
    const std::vector<AmountType>& maxSpeciesAmounts, IndexType gridSize,
    const options::IOptions& opts);

template
ReactionNetwork<NEReactionNetwork>::Cluster<plsm::OnHost>
ReactionNetwork<NEReactionNetwork>::findCluster(const Composition& comp,
    plsm::OnHost context);

template
double
ReactionNetwork<NEReactionNetwork>::getTotalConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize);

template
double
ReactionNetwork<NEReactionNetwork>::getTotalRadiusConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize);

template
double
ReactionNetwork<NEReactionNetwork>::getTotalAtomConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize);

template
double
ReactionNetwork<NEReactionNetwork>::getTotalTrappedAtomConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize);
}
}
}
