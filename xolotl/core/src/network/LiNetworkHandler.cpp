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
	NetworkType::AmountType maxV = options.getMaxV();

	// Find out how many traps we need from the options
	AmountType nTraps = 0;
	// Get the string from the options
	auto trapString = options.getTrapParameters();
	// Open the file and counts the lines if the file exist
	std::ifstream paramFile(trapString);

	if (paramFile.good()) {
		int nLines = 0;
		std::string line;

		while (std::getline(paramFile, line))
			++nLines;
		// We should have 2 lines per trap
		nTraps = nLines / 2;
	}

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		nTraps, maxH, maxV};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{nTraps + 1, maxH + 1, maxV + 1}};
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
