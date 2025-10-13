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
	NetworkType::AmountType maxV = options.getMaxV();
	NetworkType::AmountType maxI = options.getMaxI();
	NetworkType::AmountType groupingWidthV = options.getGroupingWidthA();
	// Take care of the case with no grouping
	if (options.getGroupingMin() > maxV) {
		groupingWidthV = maxV + 1;
	}
	else {
		// Adapt maxV
		int i = 0;
		while (maxV + 1 > pow(groupingWidthV, i)) {
			++i;
		}
		maxV = pow(groupingWidthV, i) - 1;
	}

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {maxV, maxI};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{groupingWidthV, maxI + 1}};
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
