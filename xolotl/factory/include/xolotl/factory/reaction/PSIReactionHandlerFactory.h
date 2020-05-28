#ifndef PSIREACTIONHANDLERFACTORY_H
#define PSIREACTIONHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/reaction/IReactionHandlerFactory.h>
#include <xolotl/core/reactants/PSIReactionNetwork.h>

namespace xolotlFactory {

/**
 * Realizes the IReactionHandlerFactory interface. Handles the network for a PSI problem.
 */
class PSIReactionHandlerFactory: public IReactionHandlerFactory {
protected:

	//! The network handler
	std::unique_ptr<xolotlCore::experimental::IReactionNetwork> theNetworkHandler;

public:

	/**
	 * The constructor creates the handlers.
	 */
	PSIReactionHandlerFactory() {
	}

	/**
	 * The destructor
	 */
	~PSIReactionHandlerFactory() {
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
		xolotlCore::experimental::PSIReactionNetwork<xolotlCore::experimental::PSIFullSpeciesList>;

		// Get the boundaries from the options
		NetworkType::AmountType maxV = options.getMaxV();
		NetworkType::AmountType maxI = options.getMaxI();
		NetworkType::AmountType maxHe =
				xolotlCore::experimental::PSIClusterGenerator<
						xolotlCore::experimental::PSIFullSpeciesList>::getMaxHePerV(
						maxV);
		NetworkType::AmountType maxD = 2.0 / 3.0 * (double) maxHe;
		NetworkType::AmountType maxT = 2.0 / 3.0 * (double) maxHe;
		NetworkType::AmountType groupingWidthHe = options.getGroupingWidthA();
		NetworkType::AmountType groupingWidthD = options.getGroupingWidthA();
		NetworkType::AmountType groupingWidthT = options.getGroupingWidthA();
		NetworkType::AmountType groupingWidthV= options.getGroupingWidthB();
		if (options.getMaxImpurity() <= 0) {
			maxHe = 0;
			groupingWidthHe = 1;
		}
		if (options.getMaxD() <= 0) {
			maxD = 0;
			groupingWidthD = 1;
		}
		if (options.getMaxT() <= 0) {
			maxT = 0;
			groupingWidthT = 1;
		}
		if (maxV <= 0) {
			maxHe = options.getMaxImpurity();
			maxD = options.getMaxD();
			maxT = options.getMaxT();
		}
		// Take care of the case with no grouping
		if (options.getGroupingMin() > maxV) {
			groupingWidthHe = 1;
			groupingWidthD = 1;
			groupingWidthT = 1;
			groupingWidthV = 1;
		}
		NetworkType::AmountType refineHe = (maxHe + 1) / groupingWidthHe;
		NetworkType::AmountType refineD = (maxD + 1) / groupingWidthD;
		NetworkType::AmountType refineT = (maxT + 1) / groupingWidthT;
		NetworkType::AmountType refineV = (maxV + 1) / groupingWidthV;
		NetworkType::AmountType refineI = (maxI + 1);

		if (maxHe + 1 != groupingWidthHe * refineHe) {
			maxHe = groupingWidthHe * (refineHe + 1) - 1;
			refineHe++;
		}
		if (maxD + 1 != groupingWidthD * refineD) {
			maxD = groupingWidthD * (refineD + 1) - 1;
			refineD++;
		}
		if (maxT + 1 != groupingWidthT * refineT) {
			maxT = groupingWidthT * (refineT + 1) - 1;
			refineT++;
		}
		if (maxV + 1 != groupingWidthV * refineV) {
			maxV = groupingWidthV * (refineV + 1) - 1;
			refineV++;
		}

		std::unique_ptr<NetworkType> rNetwork(new NetworkType( { maxHe, maxD,
				maxT, maxV, maxI }, { { refineHe, refineD, refineT, refineV,
				refineI }, { groupingWidthHe, groupingWidthD, groupingWidthT,
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

#endif // PSIREACTIONHANDLERFACTORY_H
