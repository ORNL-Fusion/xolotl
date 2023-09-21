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
	NetworkType::AmountType maxHe = options.getMaxImpurity();
	NetworkType::AmountType maxI = options.getMaxI();
	NetworkType::AmountType maxSize = options.getMaxD();
	NetworkType::AmountType groupingWidth = options.getGroupingWidthA();
	// Adapt maxSize
	int i = 0;
	while (maxSize + 1 > pow(groupingWidth, i)) {
		++i;
	}
	maxSize = pow(groupingWidth, i) - 1;
	i = 0;

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		maxHe, 1, maxSize, maxI, maxSize, maxSize, maxSize, maxI, maxSize};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{maxHe + 1, 2, groupingWidth, maxI + 1, groupingWidth, groupingWidth,
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
