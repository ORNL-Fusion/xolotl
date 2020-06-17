#include <xolotl/core/network/AlloyReactionNetwork.h>
#include <xolotl/core/network/impl/AlloyReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template
ReactionNetwork<AlloyReactionNetwork>::ReactionNetwork(
    const std::vector<AmountType>& maxSpeciesAmounts,
    const std::vector<SubdivisionRatio>& subdivisionRatios, IndexType gridSize,
    const options::IOptions& opts);

template
ReactionNetwork<AlloyReactionNetwork>::ReactionNetwork(
    const std::vector<AmountType>& maxSpeciesAmounts, IndexType gridSize,
    const options::IOptions& opts);

template
ReactionNetwork<AlloyReactionNetwork>::Cluster<plsm::OnHost>
ReactionNetwork<AlloyReactionNetwork>::findCluster(const Composition& comp,
    plsm::OnHost context);

template
double
ReactionNetwork<AlloyReactionNetwork>::getTotalConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize);

template
double
ReactionNetwork<AlloyReactionNetwork>::getTotalRadiusConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize);

template
double
ReactionNetwork<AlloyReactionNetwork>::getTotalAtomConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize);

template
double
ReactionNetwork<AlloyReactionNetwork>::getTotalTrappedAtomConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize);
}
}
}
