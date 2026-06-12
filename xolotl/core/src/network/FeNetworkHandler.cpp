#include <xolotl/core/network/FeNetworkHandler.h>
#include <xolotl/core/network/FeReactionNetwork.h>
#include <xolotl/util/Tokenizer.h>
#include <iostream> // debugging
#include <fstream>  // debugging

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
auto feNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		FeNetworkHandler>({"Fe"});
}

auto feNetworkGenerator = [](const options::IOptions& options) {
	using NetworkType = core::network::FeReactionNetwork;

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

	// Find out how many traps we need from the options
	NetworkType::AmountType nTraps = 0;
	// Get the string from the options
	auto trapString = options.getTrapParameters();
	// Break the string apart
	auto tokens = util::Tokenizer<>{trapString}();
	// We should have 4 parameters per trap
	nTraps = tokens.size() / 4;

	// Define the network
	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		nTraps, maxHe, maxV, maxI};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{nTraps + 1, groupingWidthHe, groupingWidthV, maxI + 1}};
	
	// Debugging

	std::ofstream outFile("division.logOut",std::ios::app);

	if (!outFile.is_open()) {
    		std::cerr << "Error: Could not open division.logOut\n";
	}

	// ---- Write maxSpeciesAmounts ----
	outFile << "maxSpeciesAmounts:\n";
	for (size_t i = 0; i < maxSpeciesAmounts.size(); i++) {
	    outFile << "  [" << i << "] = "
	            << maxSpeciesAmounts[i] << "\n";
	}

	// ---- Write subdivRatios ----
	outFile << "subdivRatios:\n";
	for (size_t i = 0; i < subdivRatios.size(); i++) {
	    const auto& r = subdivRatios[i];

	    outFile << "  [" << i << "] = { "
        	    << r[0] << ", "
	            << r[1] << ", "
       	 	    << r[2] << ", "
	            << r[3] << " }\n";
	}

	outFile << "-----------------------------\n";

	// Close explicitly (optional but clean)
	outFile.close();
	
	auto network = std::make_shared<NetworkType>(
                maxSpeciesAmounts, subdivRatios, 1, options);

	return network;
};

FeNetworkHandler::FeNetworkHandler(const options::IOptions& options) :
	NetworkHandler(options, feNetworkGenerator)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
