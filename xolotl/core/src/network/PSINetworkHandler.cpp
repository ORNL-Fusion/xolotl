#include <xolotl/core/Constants.h>
#include <xolotl/core/network/PSINetworkHandler.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename TSpeciesEnum,
	typename TSubdivisionRatio =
		typename PSIReactionNetwork<TSpeciesEnum>::SubdivisionRatio>
std::shared_ptr<PSIReactionNetwork<TSpeciesEnum>>
makePSIReactionNetwork(
	const std::vector<IReactionNetwork::AmountType>& maxSpeciesAmounts,
	const std::vector<TSubdivisionRatio>& subdivRatios,
	const options::IOptions& options)
{
	auto network = std::make_shared<PSIReactionNetwork<TSpeciesEnum>>(
		maxSpeciesAmounts, subdivRatios, 1, options);
	return network;
}

std::shared_ptr<IPSIReactionNetwork>
generatePSIReactionNetwork(const options::IOptions& options)
{
	using AmountType = IReactionNetwork::AmountType;

	// Get the temperature and lattice constant from the options
	double latticeConst = options.getLatticeParameter() <= 0.0 ?
		xolotl::core::tungstenLatticeConstant :
		options.getLatticeParameter();
	double temperature = options.getTempParam();

	// Get the boundaries from the options
	AmountType maxV = options.getMaxPureV();
	AmountType maxHeV = options.getMaxV();
	AmountType maxI = options.getMaxI();
	AmountType maxHe = util::getMaxHePerVLoop(maxV, latticeConst, temperature);
	AmountType maxD = 2.0 / 3.0 * (double)maxHe;
	AmountType maxT = 2.0 / 3.0 * (double)maxHe;
	AmountType groupingWidthHe = options.getGroupingWidthA();
	AmountType groupingWidthD = options.getGroupingWidthA();
	AmountType groupingWidthT = options.getGroupingWidthA();
	AmountType groupingWidthV = options.getGroupingWidthB();
	AmountType groupingWidthI = options.getGroupingWidthB();
	// To deal with pure with grouping
	AmountType deltaI = 0;

	if (options.getMaxImpurity() <= 0) {
		maxHe = 0;
		groupingWidthHe = 1;
	}
	if (options.getMaxD() <= 0) {
		maxD = 0;
		groupingWidthD = 1;
	}
	if (options.getMaxT() <= 0) {
		maxT = 0;
		groupingWidthT = 1;
	}
	if (maxV <= 0) {
		maxHe = options.getMaxImpurity();
		maxD = options.getMaxD();
		maxT = options.getMaxT();
	}
	// Take care of the case with no grouping
	if (options.getGroupingMin() > maxV) {
		groupingWidthHe = maxHe + 1;
		groupingWidthD = maxD + 1;
		groupingWidthT = maxT + 1;
		groupingWidthV = maxV + 1;
	}
	else {
		// Adapt max
		AmountType i = 0;
		while (maxHe + 1 > pow(groupingWidthHe, i)) {
			++i;
		}
		maxHe = pow(groupingWidthHe, i) - 1;
		i = 0;
		while (maxD + 1 > pow(groupingWidthD, i)) {
			++i;
		}
		maxD = pow(groupingWidthD, i) - 1;
		i = 0;
		while (maxT + 1 > pow(groupingWidthT, i)) {
			++i;
		}
		maxT = pow(groupingWidthT, i) - 1;

		// V
		i = 0;
		while (maxV + 1 > pow(groupingWidthV, i)) {
			++i;
		}
		maxV = pow(groupingWidthV, i) - 1;

		// HeV
		AmountType j = 0;
		while (maxHeV + 1 > pow(groupingWidthV, j)) {
			++j;
		}

		deltaI = i - j;
	}

	if (options.getGroupingMin() >= maxI) {
		groupingWidthI = maxI + 1;
	}
	else {
		// Adapt max
		int i = 0;
		while (maxI + 1 > pow(groupingWidthI, i)) {
			++i;
		}
		maxI = pow(groupingWidthI, i) - 1;
	}

	if (maxI > options.getGroupingMin() and maxV >= options.getGroupingMin() and
		options.getMaxImpurity() == 0) {
		// There should not be He here
		return makePSIReactionNetwork<PSIHeliumSpeciesList>(
			{0, maxV, maxI}, {{1, groupingWidthV, groupingWidthI}}, options);
	}

	if (maxD > 0 && maxT > 0) {
		std::vector<PSIReactionNetwork<PSIFullSpeciesList>::SubdivisionRatio>
			subdivRatios = {{groupingWidthHe, groupingWidthD, groupingWidthT,
				groupingWidthV, groupingWidthI}};
		auto it = subdivRatios.begin();
		for (auto k = 0; k < deltaI; k++) {
			it = subdivRatios.insert(it, {1, 1, 1, groupingWidthV, 1});
		}
		return makePSIReactionNetwork<PSIFullSpeciesList>(
			{maxHe, maxD, maxT, maxV, maxI}, subdivRatios, options);
	}
	if (maxD > 0 && maxT <= 0) {
		std::vector<
			PSIReactionNetwork<PSIDeuteriumSpeciesList>::SubdivisionRatio>
			subdivRatios = {{groupingWidthHe, groupingWidthD, groupingWidthV,
				groupingWidthI}};
		auto it = subdivRatios.begin();
		for (auto k = 0; k < deltaI; k++) {
			it = subdivRatios.insert(it, {1, 1, groupingWidthV, 1});
		}
		return makePSIReactionNetwork<PSIDeuteriumSpeciesList>(
			{maxHe, maxD, maxV, maxI}, subdivRatios, options);
	}
	if (maxD <= 0 && maxT > 0) {
		std::vector<PSIReactionNetwork<PSITritiumSpeciesList>::SubdivisionRatio>
			subdivRatios = {{groupingWidthHe, groupingWidthT, groupingWidthV,
				groupingWidthI}};
		auto it = subdivRatios.begin();
		for (auto k = 0; k < deltaI; k++) {
			it = subdivRatios.insert(it, {1, 1, groupingWidthV, 1});
		}
		return makePSIReactionNetwork<PSITritiumSpeciesList>(
			{maxHe, maxT, maxV, maxI}, subdivRatios, options);
	}
	else {
		// Either V is grouped
		if (options.getGroupingMin() > maxI) {
			AmountType refineHe = (maxHe + 1) / groupingWidthHe;
			AmountType refineV = (maxV + 1) / groupingWidthV;
			std::vector<
				PSIReactionNetwork<PSIHeliumSpeciesList>::SubdivisionRatio>
				subdivRatios = {{refineHe, refineV, maxI + 1},
					{groupingWidthHe, groupingWidthV, 1}};
			auto it = subdivRatios.begin();
			for (auto k = 0; k < deltaI; k++) {
				it = subdivRatios.insert(it, {1, groupingWidthV, 1});
			}
			return makePSIReactionNetwork<PSIHeliumSpeciesList>(
				{maxHe, maxV, maxI}, subdivRatios, options);
		}
		else {
			std::vector<
				PSIReactionNetwork<PSIHeliumSpeciesList>::SubdivisionRatio>
				subdivRatios = {
					{groupingWidthHe, groupingWidthV, groupingWidthI}};
			auto it = subdivRatios.begin();
			for (auto k = 0; k < deltaI; k++) {
				it = subdivRatios.insert(it, {1, groupingWidthV, 1});
			}
			return makePSIReactionNetwork<PSIHeliumSpeciesList>(
				{maxHe, maxV, maxI}, subdivRatios, options);
		}
	}
}

auto psiNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		PSINetworkHandler>({"W100", "W110", "W111", "W211", "Pulsed"});
} // namespace detail

PSINetworkHandler::PSINetworkHandler(const options::IOptions& options) :
	NetworkHandler(options, detail::generatePSIReactionNetwork)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
