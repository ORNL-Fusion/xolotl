#ifndef ALLOYREACTIONHANDLERFACTORY_H
#define ALLOYREACTIONHANDLERFACTORY_H

#include <memory>
#include "IReactionHandlerFactory.h"
#include <AlloyClusterNetworkLoader.h>
#include <AlloyClusterReactionNetwork.h>

namespace xolotlFactory {

/**
 * Realizes the IReactionHandlerFactory interface. Handles the network for an alloy problem.
 */
class AlloyReactionHandlerFactory: public IReactionHandlerFactory {
protected:

	//! The network loader handler
	std::shared_ptr<xolotlCore::INetworkLoader> theNetworkLoaderHandler;

	//! The network handler
	std::shared_ptr<xolotlCore::IReactionNetwork> theNetworkHandler;

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
	void initializeReactionNetwork(xolotlCore::Options &options,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
		// Get the current process ID
		int procId;
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);

		// Create a AlloyClusterNetworkLoader
		auto tempNetworkLoader = std::make_shared<
				xolotlCore::AlloyClusterNetworkLoader>(registry);
		// Give the networkFilename to the network loader
		tempNetworkLoader->setFilename(options.getNetworkFilename());
		// Set the options for the grouping scheme
		tempNetworkLoader->setMin(options.getGroupingMin());
		tempNetworkLoader->setWidth(options.getGroupingWidthA());
		theNetworkLoaderHandler = tempNetworkLoader;

		// Check if we want dummy reactions
		auto map = options.getProcesses();
		if (!map["reaction"])
			theNetworkLoaderHandler->setDummyReactions();
		// Load the network
		theNetworkHandler = theNetworkLoaderHandler->generate(options);

		if (procId == 0) {
			std::cout << "\nFactory Message: "
					<< "Master loaded network of size "
					<< theNetworkHandler->size() << "." << std::endl;
		}
	}

	/**
	 * Return the network loader.
	 *
	 * @return The network loader.
	 */
	std::shared_ptr<xolotlCore::INetworkLoader> getNetworkLoaderHandler() const {
		return theNetworkLoaderHandler;
	}

	/**
	 * Return the network.
	 *
	 * @return The network.
	 */
	std::shared_ptr<xolotlCore::IReactionNetwork> getNetworkHandler() const {
		return theNetworkHandler;
	}

};

} // end namespace xolotlFactory

#endif // ALLOYREACTIONHANDLERFACTORY_H
