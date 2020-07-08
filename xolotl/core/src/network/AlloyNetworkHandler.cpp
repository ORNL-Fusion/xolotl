#include <xolotl/core/network/AlloyNetworkHandler.h>
#include <xolotl/core/network/AlloyReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
auto alloyNetworkGenerator = [](const options::Options& options) {
	using NetworkType = AlloyReactionNetwork;

	// Get the boundaries from the options
	NetworkType::AmountType maxV = options.getMaxV();
	NetworkType::AmountType maxI = options.getMaxI();
	NetworkType::AmountType maxSize = options.getMaxImpurity();
	NetworkType::AmountType groupingWidth = options.getGroupingWidthA();
	// Adapt maxSize
	int i = 0;
	while (maxSize + 1 > pow(groupingWidth, i)) {
		++i;
	}
	maxSize = pow(groupingWidth, i) - 1;

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		maxV, maxSize, maxSize, maxI, maxSize, maxSize};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {{maxV + 1,
		groupingWidth, groupingWidth, maxI + 1, groupingWidth, groupingWidth}};
	auto network = std::make_shared<NetworkType>(
		maxSpeciesAmounts, subdivRatios, 1, options);

	network->syncClusterDataOnHost();
	network->getSubpaving().syncZones(plsm::onHost);

	return network;
};

AlloyNetworkHandler::AlloyNetworkHandler(const options::Options& options) :
	NetworkHandler(options, alloyNetworkGenerator)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
