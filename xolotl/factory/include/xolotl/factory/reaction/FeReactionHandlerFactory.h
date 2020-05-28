#ifndef FEREACTIONHANDLERFACTORY_H
#define FEREACTIONHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/reaction/IReactionHandlerFactory.h>
#include <xolotl/core/reactants/FeReactionNetwork.h>

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
			groupingWidthHe = maxHe + 1;
			groupingWidthV = maxV + 1;
		}
		else {
			// Adapt maxHe and maxV
			int i = 0;
			while (maxHe + 1 > pow(groupingWidthHe, i)) {
				++i;
			}
			maxHe = pow(groupingWidthV, i) - 1;
			i = 0;
			while (maxV + 1 > pow(groupingWidthV, i)) {
				++i;
			}
			maxV = pow(groupingWidthV, i) - 1;
		}

		std::unique_ptr<NetworkType> rNetwork(new NetworkType( { maxHe, maxV,
				maxI }, { { groupingWidthHe,
				groupingWidthV, maxI + 1 } }, 1, options));
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
