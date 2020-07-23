#include <xolotl/core/network/NENetworkHandler.h>
#include <xolotl/core/network/NEReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
auto neNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		NENetworkHandler>({"Fuel"});
}

auto neNetworkGenerator = [](const options::Options& options) {
	using NetworkType = core::network::NEReactionNetwork;

	// Get the boundaries from the options
	NetworkType::AmountType maxXe = options.getMaxImpurity();
	NetworkType::AmountType groupingWidth = options.getGroupingWidthA();
	// Take care of the case with no grouping
	if (options.getGroupingMin() > maxXe)
		groupingWidth = maxXe + 1;
	else {
		// Adapt maxXe
		int i = 0;
		while (maxXe + 1 > pow(groupingWidth, i)) {
			++i;
		}
		maxXe = pow(groupingWidth, i) - 1;
	}

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {maxXe};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {{groupingWidth}};
	// The number of grid points is set to 1 here but can be changed later
	auto network = std::make_shared<NetworkType>(
		maxSpeciesAmounts, subdivRatios, 1, options);

	network->syncClusterDataOnHost();
	network->getSubpaving().syncZones(plsm::onHost);

	return network;
};

NENetworkHandler::NENetworkHandler(const options::Options& options) :
	NetworkHandler(options, neNetworkGenerator)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
