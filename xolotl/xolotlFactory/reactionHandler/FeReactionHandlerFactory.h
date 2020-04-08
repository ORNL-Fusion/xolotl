#ifndef FEREACTIONHANDLERFACTORY_H
#define FEREACTIONHANDLERFACTORY_H

#include <memory>
#include "IReactionHandlerFactory.h"
#include <experimental/FeReactionNetwork.h>

namespace xolotlFactory {

/**
 * Realizes the IReactionHandlerFactory interface. Handles the network for an iron problem.
 */
class FeReactionHandlerFactory: public IReactionHandlerFactory {
protected:

	//! The network handler
	std::unique_ptr<xolotlCore::experimental::IReactionNetwork> theNetworkHandler;

public:

	/**
	 * The constructor creates the handlers.
	 */
	FeReactionHandlerFactory() {
	}

	/**
	 * The destructor
	 */
	~FeReactionHandlerFactory() {
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
		xolotlCore::experimental::FeReactionNetwork;

		// Get the boundaries from the options
		NetworkType::AmountType maxV = options.getMaxV();
		NetworkType::AmountType maxI = options.getMaxI();
		NetworkType::AmountType maxHe = options.getMaxImpurity();
		NetworkType::AmountType groupingWidthHe = options.getGroupingWidthA();
		NetworkType::AmountType groupingWidthV = options.getGroupingWidthB();
		// Take care of the case with no grouping
		if (options.getGroupingMin() > maxV) {
			groupingWidthHe = 1;
			groupingWidthV = 1;
		}
		NetworkType::AmountType refineHe = (maxHe + 1) / groupingWidthHe;
		NetworkType::AmountType refineV = (maxV + 1) / groupingWidthV;
		NetworkType::AmountType refineI = (maxI + 1);

		if (maxHe + 1 != groupingWidthHe * refineHe) {
			maxHe = groupingWidthHe * (refineHe + 1) - 1;
			refineHe++;
		}
		if (maxV + 1 != groupingWidthV * refineV) {
			maxV = groupingWidthV * (refineV + 1) - 1;
			refineV++;
		}

		std::unique_ptr<NetworkType> rNetwork(new NetworkType( { maxHe, maxV,
				maxI }, { { refineHe, refineV, refineI }, { groupingWidthHe,
				groupingWidthV, 1 } }, 1, options));
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

#endif // FEREACTIONHANDLERFACTORY_H
