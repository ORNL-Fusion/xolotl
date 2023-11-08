#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/core/network/impl/PSIReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template class PSIReactionNetwork<PSIFullSpeciesList>;

template ReactionNetwork<PSIReactionNetwork<PSIFullSpeciesList>>::
	ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
		const std::vector<SubdivisionRatio>& subdivisionRatios,
		IndexType gridSize, const options::IOptions& opts);

template ReactionNetwork<PSIReactionNetwork<PSIFullSpeciesList>>::
	ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
		IndexType gridSize, const options::IOptions& opts);

template double
ReactionNetwork<PSIReactionNetwork<PSIFullSpeciesList>>::getTotalConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIFullSpeciesList>>::
	getTotalRadiusConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIFullSpeciesList>>::
	getTotalAtomConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIFullSpeciesList>>::
	getTotalTrappedAtomConcentration(
		ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIFullSpeciesList>>::getTotalVolumeFraction(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<PSIReactionNetwork<PSIFullSpeciesList>>::getSmallConcentration(
	ConcentrationsView concentrations, Species type, AmountType maxSize);
} // namespace network
} // namespace core
} // namespace xolotl
