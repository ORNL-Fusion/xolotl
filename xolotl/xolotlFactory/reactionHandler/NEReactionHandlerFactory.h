#ifndef NEREACTIONHANDLERFACTORY_H
#define NEREACTIONHANDLERFACTORY_H

#include <memory>
#include "IReactionHandlerFactory.h"
#include <experimental/NEReactionNetwork.h>

namespace xolotlFactory {

/**
 * Realizes the IReactionHandlerFactory interface. Handles the network for a NE problem.
 */
class NEReactionHandlerFactory: public IReactionHandlerFactory {
protected:

	//! The network handler
	std::unique_ptr<xolotlCore::experimental::IReactionNetwork> theNetworkHandler;

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

		using NetworkType =
		xolotlCore::experimental::NEReactionNetwork;

		// Get the boundaries from the options
		NetworkType::AmountType maxXe = options.getMaxImpurity();
		NetworkType::AmountType groupingWidth = options.getGroupingWidthA();
		NetworkType::AmountType refine = (maxXe + 1) / groupingWidth;


		std::unique_ptr<NetworkType> rNetwork(
					new NetworkType( { maxXe }, { { refine }, { groupingWidth } }, 1,
				options));
		rNetwork->syncClusterDataOnHost();
		rNetwork->getSubpaving().syncZones(plsm::onHost);
		theNetworkHandler = std::move( rNetwork );

		if (procId == 0) {
			std::cout << "\nFactory Message: "
					<< "Master loaded network of size "
					<< theNetworkHandler->getDOF() << "." << std::endl;
		}
//		// Set the fission rate in the network to compute the diffusion coefficient correctly
//		theNetworkHandler->setFissionRate(options.getFluxAmplitude());
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

#endif // NEREACTIONHANDLERFACTORY_H
