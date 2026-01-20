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

auto neNetworkGenerator = [](const options::IOptions& options) {
	using NetworkType = core::network::NEReactionNetwork;

	// Get the boundaries from the options
	NetworkType::AmountType maxV = options.getMaxV();
	NetworkType::AmountType maxI = options.getMaxI();
	NetworkType::AmountType maxXe = options.getMaxImpurity();
	NetworkType::AmountType groupingWidthXe = options.getGroupingWidthA();
	NetworkType::AmountType groupingWidthV = options.getGroupingWidthB();
	// Take care of the case with no grouping
	if (options.getGroupingMin() > maxV) {
		groupingWidthXe = maxXe + 1;
		groupingWidthV = maxV + 1;
	}
	else {
		// Adapt maxXe and maxV
		int i = 0;
		while (maxXe + 1 > pow(groupingWidthXe, i)) {
			++i;
		}
		maxXe = pow(groupingWidthXe, i) - 1;
		i = 0;
		while (maxV + 1 > pow(groupingWidthV, i)) {
			++i;
		}
		maxV = pow(groupingWidthV, i) - 1;
	}

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		maxXe, maxV, maxI};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{groupingWidthXe, groupingWidthV, maxI + 1}};
	auto network = std::make_shared<NetworkType>(
		maxSpeciesAmounts, subdivRatios, 1, options);

	return network;
};

NENetworkHandler::NENetworkHandler(const options::IOptions& options) :
	NetworkHandler(options, neNetworkGenerator)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
