#ifndef NEREACTIONHANDLERFACTORY_H
#define NEREACTIONHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/network/IReactionHandlerFactory.h>
#include <xolotl/core/network/NEReactionNetwork.h>

namespace xolotl {
namespace factory {
namespace network {

/**
 * Realizes the IReactionHandlerFactory interface. Handles the network for a NE problem.
 */
class NEReactionHandlerFactory: public IReactionHandlerFactory {
protected:

	//! The network handler
	std::unique_ptr<core::network::IReactionNetwork> theNetworkHandler;

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
	void initializeReactionNetwork(const options::Options &opts,
			std::shared_ptr<perf::IHandlerRegistry> registry) {
		// Get the current process ID
		int procId;
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);

		using NetworkType =
		core::network::NEReactionNetwork;

		// Get the boundaries from the options
		NetworkType::AmountType maxXe = opts.getMaxImpurity();
		NetworkType::AmountType groupingWidth = opts.getGroupingWidthA();
		// Take care of the case with no grouping
		if (opts.getGroupingMin() > maxXe)
			groupingWidth = maxXe + 1;
		else {
			// Adapt maxXe
			int i = 0;
			while (maxXe + 1 > pow(groupingWidth, i)) {
				++i;
			}
			maxXe = pow(groupingWidth, i) - 1;
		}

		// The number of grid points is set to 1 here but can be changed later
		std::unique_ptr<NetworkType> rNetwork(new NetworkType( { maxXe }, { {
				groupingWidth } }, 1, opts));
		rNetwork->syncClusterDataOnHost();
		rNetwork->getSubpaving().syncZones(plsm::onHost);
		theNetworkHandler = std::move(rNetwork);

		if (procId == 0) {
			std::cout << "\nFactory Message: "
					<< "Master loaded network of size "
					<< theNetworkHandler->getDOF() << "." << std::endl;
		}
//		// Set the fission rate in the network to compute the diffusion coefficient correctly
//		theNetworkHandler->setFissionRate(opts.getFluxAmplitude());
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

#endif // NEREACTIONHANDLERFACTORY_H
