#include <xolotl/core/network/AlloyNetworkHandler.h>
#include <xolotl/core/network/AlloyReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
auto alloyNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		AlloyNetworkHandler>({"800H"});
}

auto alloyNetworkGenerator = [](const options::IOptions& options) {
	using NetworkType = AlloyReactionNetwork;

	// Get the boundaries from the options
	NetworkType::AmountType maxV = options.getMaxV();
	NetworkType::AmountType maxI = options.getMaxI();
	NetworkType::AmountType maxSize = options.getMaxImpurity();
	NetworkType::AmountType maxVoid = options.getMaxD();
	NetworkType::AmountType groupingWidth = options.getGroupingWidthA();
	// Adapt maxSize
	int i = 0;
	while (maxSize + 1 > pow(groupingWidth, i)) {
		++i;
	}
	maxSize = pow(groupingWidth, i) - 1;
	i = 0;
	while (maxVoid + 1 > pow(groupingWidth, i)) {
		++i;
	}
	maxVoid = pow(groupingWidth, i) - 1;

	NetworkType::AmountType maxVLoopSize = maxV == 0 ? 0 : maxSize;
	NetworkType::AmountType maxILoopSize = maxI == 0 ? 0 : maxSize;
	NetworkType::AmountType groupingV = maxV == 0 ? 1 : groupingWidth;
	NetworkType::AmountType groupingI = maxI == 0 ? 1 : groupingWidth;

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		maxV, maxVoid, maxVLoopSize, maxI, maxILoopSize, maxILoopSize};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{maxV + 1, groupingV, groupingV, maxI + 1, groupingI, groupingI}};
	auto network = std::make_shared<NetworkType>(
		maxSpeciesAmounts, subdivRatios, 1, options);

	return network;
};

AlloyNetworkHandler::AlloyNetworkHandler(const options::IOptions& options) :
	NetworkHandler(options, alloyNetworkGenerator)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
