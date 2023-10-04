#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/core/network/impl/PSIReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template class PSIReactionNetwork<PSIHeliumSpeciesList>;

template ReactionNetwork<PSIReactionNetwork<PSIHeliumSpeciesList>>::
	ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
		const std::vector<SubdivisionRatio>& subdivisionRatios,
		IndexType gridSize, const options::IOptions& opts);

template ReactionNetwork<PSIReactionNetwork<PSIHeliumSpeciesList>>::
	ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
		IndexType gridSize, const options::IOptions& opts);

template double
ReactionNetwork<PSIReactionNetwork<PSIHeliumSpeciesList>>::
	getTotalConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIHeliumSpeciesList>>::
	getTotalRadiusConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIHeliumSpeciesList>>::
	getTotalAtomConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIHeliumSpeciesList>>::
	getTotalTrappedAtomConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIHeliumSpeciesList>>::
	getTotalVolumeFraction(
		ConcentrationsView concentrations, Species type, AmountType minSize);
} // namespace network
} // namespace core
} // namespace xolotl
