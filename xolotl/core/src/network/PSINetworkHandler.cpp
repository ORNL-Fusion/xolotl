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
	network->syncClusterDataOnHost();
	network->getSubpaving().syncZones(plsm::onHost);
	return network;
}

std::shared_ptr<IPSIReactionNetwork>
generatePSIReactionNetwork(const options::Options& options)
{
	using AmountType = IReactionNetwork::AmountType;

	// Get the boundaries from the options
	AmountType maxV = options.getMaxV();
	AmountType maxI = options.getMaxI();
	AmountType maxHe = psi::getMaxHePerV(maxV, options.getHeVRatio());
	AmountType maxD = 2.0 / 3.0 * (double)maxHe;
	AmountType maxT = 2.0 / 3.0 * (double)maxHe;
	AmountType groupingWidthHe = options.getGroupingWidthA();
	AmountType groupingWidthD = options.getGroupingWidthA();
	AmountType groupingWidthT = options.getGroupingWidthA();
	AmountType groupingWidthV = options.getGroupingWidthB();
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
		groupingWidthHe = 1;
		groupingWidthD = 1;
		groupingWidthT = 1;
		groupingWidthV = 1;
	}
	AmountType refineHe = (maxHe + 1) / groupingWidthHe;
	AmountType refineD = (maxD + 1) / groupingWidthD;
	AmountType refineT = (maxT + 1) / groupingWidthT;
	AmountType refineV = (maxV + 1) / groupingWidthV;
	AmountType refineI = (maxI + 1);

	if (maxHe + 1 != groupingWidthHe * refineHe) {
		maxHe = groupingWidthHe * (refineHe + 1) - 1;
		refineHe++;
	}
	if (maxD > 0) {
		if (maxD + 1 != groupingWidthD * refineD) {
			maxD = groupingWidthD * (refineD + 1) - 1;
			refineD++;
		}
	}
	if (maxT > 0) {
		if (maxT + 1 != groupingWidthT * refineT) {
			maxT = groupingWidthT * (refineT + 1) - 1;
			refineT++;
		}
	}
	if (maxV + 1 != groupingWidthV * refineV) {
		maxV = groupingWidthV * (refineV + 1) - 1;
		refineV++;
	}

	if (maxD > 0 && maxT > 0) {
		return makePSIReactionNetwork<PSIFullSpeciesList>(
			{maxHe, maxD, maxT, maxV, maxI},
			{{refineHe, refineD, refineT, refineV, refineI},
				{groupingWidthHe, groupingWidthD, groupingWidthT,
					groupingWidthV, 1}},
			options);
	}
	else if (maxD > 0 && maxT <= 0) {
        return makePSIReactionNetwork<PSIDeuteriumSpeciesList>(
            {maxHe, maxD, maxV, maxI},
            {{refineHe, refineD, refineV, refineI},
                {groupingWidthHe, groupingWidthD, groupingWidthV, 1}},
            options);
	}
	else if (maxD <= 0 && maxT > 0) {
        return makePSIReactionNetwork<PSITritiumSpeciesList>(
            {maxHe, maxT, maxV, maxI},
            {{refineHe, refineT, refineV, refineI},
                {groupingWidthHe, groupingWidthT, groupingWidthV, 1}},
            options);
	}
	else {
		return makePSIReactionNetwork<PSIHeliumSpeciesList>({maxHe, maxV, maxI},
			{{refineHe, refineV, refineI},
				{groupingWidthHe, groupingWidthV, 1}},
			options);
	}
}

auto psiNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		PSINetworkHandler>({"W100", "W110", "W111", "W211", "Pulsed"});
} // namespace detail

PSINetworkHandler::PSINetworkHandler(const options::Options& options) :
	NetworkHandler(options, detail::generatePSIReactionNetwork)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
