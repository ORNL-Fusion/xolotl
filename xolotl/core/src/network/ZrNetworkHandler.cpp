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
	// std::cout << "NetworkHandler-start \n";

	using NetworkType = ZrReactionNetwork;

	// Get the boundaries from the options
	NetworkType::AmountType maxV = options.getMaxV();
	NetworkType::AmountType maxI = options.getMaxI();

	int i = 0;
	while (maxV + 1 > pow(2, i)) {
		++i;
	}
	maxV = pow(2, i) - 1;
	i = 0;
	while (maxI + 1 > pow(2, i)) {
		++i;
	}
	maxI = pow(2, i) - 1;

	// adding basal
	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {maxV, maxV, maxI};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {{2, 2, 2}};

	auto network = std::make_shared<NetworkType>(
		maxSpeciesAmounts, subdivRatios, 1, options);

	network->syncClusterDataOnHost();

	return network;
};

ZrNetworkHandler::ZrNetworkHandler(const options::IOptions& options) :
	NetworkHandler(options, zrNetworkGenerator)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
