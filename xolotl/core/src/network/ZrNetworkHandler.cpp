#include <xolotl/core/network/ZrNetworkHandler.h>
#include <xolotl/core/network/ZrReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
auto zrNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		ZrNetworkHandler>({"AlphaZr"});
}

auto zrNetworkGenerator = [](const options::IOptions& options) {
	using NetworkType = ZrReactionNetwork;

	// Get the boundaries from the options
	NetworkType::AmountType maxV = options.getMaxV();
	NetworkType::AmountType maxI = options.getMaxI();

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {maxV, maxI};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{maxV + 1, maxI + 1}};
	auto network = std::make_shared<NetworkType>(
		maxSpeciesAmounts, subdivRatios, 1, options);

	network->syncClusterDataOnHost();
	network->getSubpaving().syncZones(plsm::onHost);

	return network;
};

ZrNetworkHandler::ZrNetworkHandler(const options::IOptions& options) :
	NetworkHandler(options, zrNetworkGenerator)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
