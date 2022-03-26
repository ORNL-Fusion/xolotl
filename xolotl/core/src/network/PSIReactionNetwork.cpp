#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/core/network/impl/PSIReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
// PSIFullSpeciesList
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

// PSIHeliumSpeciesList
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

// PSIDeuteriumSpeciesList
template class PSIReactionNetwork<PSIDeuteriumSpeciesList>;

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

// PSITritiumSpeciesList
template class PSIReactionNetwork<PSITritiumSpeciesList>;

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
