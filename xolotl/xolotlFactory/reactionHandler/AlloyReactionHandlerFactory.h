#ifndef ALLOYREACTIONHANDLERFACTORY_H
#define ALLOYREACTIONHANDLERFACTORY_H

#include <memory>
#include "IReactionHandlerFactory.h"
#include <experimental/AlloyReactionNetwork.h>

namespace xolotlFactory {

/**
 * Realizes the IReactionHandlerFactory interface. Handles the network for an alloy problem.
 */
class AlloyReactionHandlerFactory: public IReactionHandlerFactory {
protected:

	//! The network handler
	std::unique_ptr<xolotlCore::experimental::IReactionNetwork> theNetworkHandler;

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
	 * @param options The options.
	 * @param registry The performance registry.
	 */
	void initializeReactionNetwork(const xolotlCore::Options &options,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
		// Get the current process ID
		int procId;
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);

		using NetworkType =
		xolotlCore::experimental::AlloyReactionNetwork;

		// Get the boundaries from the options
		NetworkType::AmountType maxV = options.getMaxV();
		NetworkType::AmountType maxI = options.getMaxI();
		NetworkType::AmountType maxSize = options.getMaxImpurity();
		NetworkType::AmountType groupingWidth = options.getGroupingWidthA();
		// Take care of the case with no grouping
		if (options.getGroupingMin() > maxSize) {
			groupingWidth = maxSize + 1;
		}

		std::unique_ptr<NetworkType> rNetwork(
				new NetworkType( { maxV, maxSize, maxSize, maxI, maxSize,
						maxSize }, { { maxV + 1, groupingWidth,
						groupingWidth, maxI + 1, groupingWidth, groupingWidth } }, 1,
						options));
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
	xolotlCore::experimental::IReactionNetwork& getNetworkHandler() const {
		return *theNetworkHandler;
	}

};

} // end namespace xolotlFactory

#endif // ALLOYREACTIONHANDLERFACTORY_H
