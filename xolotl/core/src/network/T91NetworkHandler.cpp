#include <xolotl/core/network/T91NetworkHandler.h>
#include <xolotl/core/network/T91ReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
auto t91NetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		T91NetworkHandler>({"T91"});
}

auto t91NetworkGenerator = [](const options::IOptions& options) {
	using NetworkType = core::network::T91ReactionNetwork;

	// Get the boundaries from the options
	NetworkType::AmountType maxV = options.getMaxV();
	NetworkType::AmountType maxI = options.getMaxI();
	NetworkType::AmountType maxHe = options.getMaxImpurity();
	NetworkType::AmountType groupingWidthHe = options.getGroupingWidthA();
	NetworkType::AmountType groupingWidthV = options.getGroupingWidthB();
	// Take care of the case with no grouping
	if (options.getGroupingMin() > maxV) {
		groupingWidthHe = maxHe + 1;
		groupingWidthV = maxV + 1;
	}
	else {
		// Adapt maxHe and maxV
		int i = 0;
		while (maxHe + 1 > pow(groupingWidthHe, i)) {
			++i;
		}
		maxHe = pow(groupingWidthHe, i) - 1;
		i = 0;
		while (maxV + 1 > pow(groupingWidthV, i)) {
			++i;
		}
		maxV = pow(groupingWidthV, i) - 1;
	}

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		maxHe, maxV, maxI};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{groupingWidthHe, groupingWidthV, maxI + 1}};
	auto network = std::make_shared<NetworkType>(
		maxSpeciesAmounts, subdivRatios, 1, options);

	return network;
};

T91NetworkHandler::T91NetworkHandler(const options::IOptions& options) :
	NetworkHandler(options, t91NetworkGenerator)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
