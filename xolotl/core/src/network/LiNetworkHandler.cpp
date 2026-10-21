#include <xolotl/core/network/LiNetworkHandler.h>
#include <xolotl/core/network/LiReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
auto liNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		LiNetworkHandler>({"Li"});
}

auto liNetworkGenerator = [](const options::IOptions& options) {
	using NetworkType = core::network::LiReactionNetwork;

	// Get the boundaries from the options
	NetworkType::AmountType maxH = options.getMaxImpurity();

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {maxH};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {{maxH + 1}};
	auto network = std::make_shared<NetworkType>(
		maxSpeciesAmounts, subdivRatios, 1, options);

	return network;
};

LiNetworkHandler::LiNetworkHandler(const options::IOptions& options) :
	NetworkHandler(options, liNetworkGenerator)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
