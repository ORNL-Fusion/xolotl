#include <xolotl/core/network/RPVNetworkHandler.h>
#include <xolotl/core/network/RPVReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
auto rpvNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		RPVNetworkHandler>({"RPV"});
}

auto rpvNetworkGenerator = [](const options::IOptions& options) {
	using NetworkType = core::network::RPVReactionNetwork;

	// Get the boundaries from the options
	NetworkType::AmountType maxSize = options.getMaxV();
	NetworkType::AmountType groupingWidth = options.getGroupingWidthA();
	// Take care of the case with no grouping
	if (options.getGroupingMin() > maxSize) {
		groupingWidth = maxSize + 1;
	}
	else {
		// Adapt maxSize
		int i = 0;
		while (maxSize + 1 > pow(groupingWidth, i)) {
			++i;
		}
		maxSize = pow(groupingWidth, i) - 1;
	}

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		maxSize, maxSize, maxSize, maxSize};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{groupingWidth, groupingWidth, groupingWidth, groupingWidth}};
	auto network = std::make_shared<NetworkType>(
		maxSpeciesAmounts, subdivRatios, 1, options);

	return network;
};

RPVNetworkHandler::RPVNetworkHandler(const options::IOptions& options) :
	NetworkHandler(options, rpvNetworkGenerator)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
