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
	NetworkType::AmountType maxV = options.getMaxImpurity();
	NetworkType::AmountType groupingWidth = options.getGroupingWidthA();
	// Adapt maxSize
	int i = 0;
	while (maxV + 1 > pow(groupingWidth, i)) {
		++i;
	}
	maxV = pow(groupingWidth, i) - 1;
	i = 0;
	while (maxI + 1 > pow(groupingWidth, i)) {
		++i;
	}
	maxI = pow(groupingWidth, i) - 1;

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		0, maxV, maxI, 0, 0, 0, 0, 0};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{1, groupingWidth, groupingWidth, 1, 1, 1, 1, 1}};
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
