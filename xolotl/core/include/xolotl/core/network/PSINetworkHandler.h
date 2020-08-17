#pragma once

#include <xolotl/core/network/NetworkHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename TSpeciesEnum>
auto
generatePSIReactionNetwork(const options::Options& options)
{
	using NetworkType = PSIReactionNetwork<TSpeciesEnum>;
	using AmountType = typename NetworkType::AmountType;

	// Get the boundaries from the options
	AmountType maxV = options.getMaxV();
	AmountType maxI = options.getMaxI();
	AmountType maxHe = PSIClusterGenerator<TSpeciesEnum>::getMaxHePerV(
		maxV, options.getHeVRatio());
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
	if (maxD + 1 != groupingWidthD * refineD) {
		maxD = groupingWidthD * (refineD + 1) - 1;
		refineD++;
	}
	if (maxT + 1 != groupingWidthT * refineT) {
		maxT = groupingWidthT * (refineT + 1) - 1;
		refineT++;
	}
	if (maxV + 1 != groupingWidthV * refineV) {
		maxV = groupingWidthV * (refineV + 1) - 1;
		refineV++;
	}

	std::vector<AmountType> maxSpeciesAmounts = {maxHe, maxD, maxT, maxV, maxI};
	std::vector<typename NetworkType::SubdivisionRatio> subdivRatios = {
		{refineHe, refineD, refineT, refineV, refineI},
		{groupingWidthHe, groupingWidthD, groupingWidthT, groupingWidthV, 1}};
	auto network = std::make_shared<NetworkType>(
		maxSpeciesAmounts, subdivRatios, 1, options);

	network->syncClusterDataOnHost();
	network->getSubpaving().syncZones(plsm::onHost);

	return network;
}
} // namespace detail

template <typename TSpeciesEnum>
class PSINetworkHandler : public NetworkHandler
{
public:
	PSINetworkHandler(const options::Options& options) :
		NetworkHandler(
			options, detail::generatePSIReactionNetwork<TSpeciesEnum>)
	{
	}
};
} // namespace network
} // namespace core
} // namespace xolotl
