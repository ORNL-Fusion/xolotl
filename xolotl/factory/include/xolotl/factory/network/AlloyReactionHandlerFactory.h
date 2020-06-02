#ifndef ALLOYREACTIONHANDLERFACTORY_H
#define ALLOYREACTIONHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/network/IReactionHandlerFactory.h>
#include <xolotl/core/network/AlloyReactionNetwork.h>

namespace xolotl {
namespace factory {
namespace network {

/**
 * Realizes the IReactionHandlerFactory interface. Handles the network for an alloy problem.
 */
class AlloyReactionHandlerFactory: public IReactionHandlerFactory {
protected:

	//! The network handler
	std::unique_ptr<core::network::IReactionNetwork> theNetworkHandler;

public:

	/**
	 * The constructor creates the handlers.
	 */
	AlloyReactionHandlerFactory() {
	}

	/**
	 * The destructor
	 */
	~AlloyReactionHandlerFactory() {
	}

	/**
	 * Initialize the reaction network.
	 *
	 * @param opts The options.
	 * @param registry The performance registry.
	 */
	void initializeReactionNetwork(const options::Options &opts,
			std::shared_ptr<perf::IHandlerRegistry> registry) {
		// Get the current process ID
		int procId;
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);

		using NetworkType =
		core::network::AlloyReactionNetwork;

		// Get the boundaries from the options
		NetworkType::AmountType maxV = opts.getMaxV();
		NetworkType::AmountType maxI = opts.getMaxI();
		NetworkType::AmountType maxSize = opts.getMaxImpurity();
		NetworkType::AmountType groupingWidth = opts.getGroupingWidthA();
		// Take care of the case with no grouping
		if (opts.getGroupingMin() > maxSize) {
			groupingWidth = maxSize + 1;
		}

		std::unique_ptr<NetworkType> rNetwork(
				new NetworkType( { maxV, maxSize, maxSize, maxI, maxSize,
						maxSize }, { { maxV + 1, groupingWidth,
						groupingWidth, maxI + 1, groupingWidth, groupingWidth } }, 1,
						opts));
		rNetwork->syncClusterDataOnHost();
		rNetwork->getSubpaving().syncZones(plsm::onHost);
		theNetworkHandler = std::move(rNetwork);

		if (procId == 0) {
			std::cout << "\nFactory Message: "
					<< "Master loaded network of size "
					<< theNetworkHandler->getDOF() << "." << std::endl;
		}
	}

	/**
	 * Return the network.
	 *
	 * @return The network.
	 */
    core::network::IReactionNetwork& getNetworkHandler() const {
		return *theNetworkHandler;
	}

};

} // end namespace network
} // end namespace factory
} // end namespace xolotl

#endif // ALLOYREACTIONHANDLERFACTORY_H
