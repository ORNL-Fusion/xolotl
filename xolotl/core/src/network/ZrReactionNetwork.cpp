#include <xolotl/core/network/ZrReactionNetwork.h>
#include <xolotl/core/network/impl/ZrReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template ReactionNetwork<ZrReactionNetwork>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts,
	const std::vector<SubdivisionRatio>& subdivisionRatios, IndexType gridSize,
	const options::IOptions& opts);

template ReactionNetwork<ZrReactionNetwork>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts, IndexType gridSize,
	const options::IOptions& opts);

template double
ReactionNetwork<ZrReactionNetwork>::getTotalConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<ZrReactionNetwork>::getTotalRadiusConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<ZrReactionNetwork>::getTotalAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<ZrReactionNetwork>::getTotalTrappedAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<ZrReactionNetwork>::getTotalVolumeFraction(
	ConcentrationsView concentrations, Species type, AmountType minSize);
} // namespace network
} // namespace core
} // namespace xolotl
