#include <xolotl/core/network/FeCrNetworkHandler.h>
#include <xolotl/core/network/FeCrReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
auto fecrNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		FeCrNetworkHandler>({"FeCr"});
}

auto fecrNetworkGenerator = [](const options::IOptions& options) {
	using NetworkType = FeCrReactionNetwork;

	// Get the boundaries from the options
	NetworkType::AmountType maxI = options.getMaxI();
	NetworkType::AmountType maxSize = options.getMaxImpurity();
	//	NetworkType::AmountType groupingWidth = options.getGroupingWidthA();
	NetworkType::AmountType groupingWidth = 2;
	// Adapt maxSize
	int i = 0;
	while (maxSize + 1 > pow(groupingWidth, i)) {
		++i;
	}
	maxSize = pow(groupingWidth, i) - 1;

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		1, maxSize, maxI, maxSize, maxSize, maxSize, maxI, maxSize};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{2, groupingWidth, maxI + 1, groupingWidth, groupingWidth,
			groupingWidth, maxI + 1, groupingWidth}};
	auto network = std::make_shared<NetworkType>(
		maxSpeciesAmounts, subdivRatios, 1, options);

	network->syncClusterDataOnHost();

	return network;
};

FeCrNetworkHandler::FeCrNetworkHandler(const options::IOptions& options) :
	NetworkHandler(options, fecrNetworkGenerator)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
