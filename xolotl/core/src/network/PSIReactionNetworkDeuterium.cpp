#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/core/network/impl/PSIReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template class PSIReactionNetwork<PSIDeuteriumSpeciesList>;

template ReactionNetwork<
	PSIReactionNetwork<PSIDeuteriumSpeciesList>>::ReactionNetwork();

template ReactionNetwork<PSIReactionNetwork<PSIDeuteriumSpeciesList>>::
	ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
		const std::vector<SubdivisionRatio>& subdivisionRatios,
		IndexType gridSize, const options::IOptions& opts);

template ReactionNetwork<PSIReactionNetwork<PSIDeuteriumSpeciesList>>::
	ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
		IndexType gridSize, const options::IOptions& opts);

template double
ReactionNetwork<PSIReactionNetwork<PSIDeuteriumSpeciesList>>::
	getTotalConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIDeuteriumSpeciesList>>::
	getTotalRadiusConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIDeuteriumSpeciesList>>::
	getTotalAtomConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIDeuteriumSpeciesList>>::
	getTotalTrappedAtomConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIDeuteriumSpeciesList>>::
	getTotalVolumeFraction(
		ConcentrationsView concentrations, Species type, AmountType minSize);
} // namespace network
} // namespace core
} // namespace xolotl
