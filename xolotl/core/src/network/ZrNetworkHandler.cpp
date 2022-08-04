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
	NetworkType::AmountType maxB = options.getMaxImpurity();
	NetworkType::AmountType maxV = options.getMaxV();
	NetworkType::AmountType maxI = options.getMaxI();
	NetworkType::AmountType ratioB = 1;
	NetworkType::AmountType ratioV = 1;
	NetworkType::AmountType ratioI = 1;

	if (maxV > 0) {
		int i = 0;
		while (maxV + 1 > pow(2, i)) {
			++i;
		}
		maxV = pow(2, i) - 1;
		ratioV = 2;
	}
	if (maxI > 0) {
		int i = 0;
		while (maxI + 1 > pow(2, i)) {
			++i;
		}
		maxI = pow(2, i) - 1;
		ratioI = 2;
	}
	if (maxB > 0) {
		int i = 0;
		while (maxB + 1 > pow(2, i)) {
			++i;
		}
		maxB = pow(2, i) - 1;
		ratioB = 2;
	}

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
	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {maxV, maxB, maxI};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{ratioV, ratioB, ratioI}};

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
