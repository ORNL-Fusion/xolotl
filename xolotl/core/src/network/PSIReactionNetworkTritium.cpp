#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/core/network/impl/PSIReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template class PSIReactionNetwork<PSITritiumSpeciesList>;

template ReactionNetwork<
	PSIReactionNetwork<PSITritiumSpeciesList>>::ReactionNetwork();

template ReactionNetwork<PSIReactionNetwork<PSITritiumSpeciesList>>::
	ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
		const std::vector<SubdivisionRatio>& subdivisionRatios,
		IndexType gridSize, const options::IOptions& opts);

template ReactionNetwork<PSIReactionNetwork<PSITritiumSpeciesList>>::
	ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
		IndexType gridSize, const options::IOptions& opts);

template double
ReactionNetwork<PSIReactionNetwork<PSITritiumSpeciesList>>::
	getTotalConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSITritiumSpeciesList>>::
	getTotalRadiusConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSITritiumSpeciesList>>::
	getTotalAtomConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSITritiumSpeciesList>>::
	getTotalTrappedAtomConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSITritiumSpeciesList>>::
	getTotalVolumeFraction(
		ConcentrationsView concentrations, Species type, AmountType minSize);
} // namespace network
} // namespace core
} // namespace xolotl
