#ifndef NEREACTIONHANDLERFACTORY_H
#define NEREACTIONHANDLERFACTORY_H

#include <memory>
#include "IReactionHandlerFactory.h"
#include <NEClusterNetworkLoader.h>
#include <NEClusterReactionNetwork.h>

namespace xolotlFactory {

/**
 * Realizes the IReactionHandlerFactory interface. Handles the network for a NE problem.
 */
class NEReactionHandlerFactory: public IReactionHandlerFactory {
protected:

	//! The network loader handler
	std::shared_ptr<xolotlCore::INetworkLoader> theNetworkLoaderHandler;

	//! The network handler
	std::unique_ptr<xolotlCore::IReactionNetwork> theNetworkHandler;

public:

	/**
	 * The constructor creates the handlers.
	 */
	NEReactionHandlerFactory() {
	}

	/**
	 * The destructor
	 */
	~NEReactionHandlerFactory() {
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

		// Create a NEClusterNetworkLoader
		auto tempNetworkLoader = std::make_shared<
				xolotlCore::NEClusterNetworkLoader>(registry);
		// Give the networkFilename to the network loader
		tempNetworkLoader->setFilename(options.getNetworkFilename());
		// Set the options for the grouping scheme
		tempNetworkLoader->setXeMin(options.getGroupingMin());
		tempNetworkLoader->setWidth(options.getGroupingWidthA());
		theNetworkLoaderHandler = tempNetworkLoader;

		// Check if we want dummy reactions
		auto map = options.getProcesses();
		if (!map["reaction"])
			theNetworkLoaderHandler->setDummyReactions();
		// Load the network
		if (options.useHDF5())
			theNetworkHandler = theNetworkLoaderHandler->load(options);
		else
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
	xolotlCore::IReactionNetwork& getNetworkHandler() const {
		return *theNetworkHandler;
	}

};

} // end namespace xolotlFactory

#endif // NEREACTIONHANDLERFACTORY_H
