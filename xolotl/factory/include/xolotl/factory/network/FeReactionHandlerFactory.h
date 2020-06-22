#ifndef FEREACTIONHANDLERFACTORY_H
#define FEREACTIONHANDLERFACTORY_H

#include <memory>

#include <xolotl/core/network/FeReactionNetwork.h>
#include <xolotl/factory/network/IReactionHandlerFactory.h>

namespace xolotl
{
namespace factory
{
namespace network
{
/**
 * Realizes the IReactionHandlerFactory interface. Handles the network for an
 * iron problem.
 */
class FeReactionHandlerFactory : public IReactionHandlerFactory
{
protected:
	//! The network handler
	std::unique_ptr<core::network::IReactionNetwork> theNetworkHandler;

public:
	/**
	 * The constructor creates the handlers.
	 */
	FeReactionHandlerFactory()
	{
	}

	/**
	 * The destructor
	 */
	~FeReactionHandlerFactory()
	{
	}

	/**
	 * Initialize the reaction network.
	 *
	 * @param opts The options.
	 * @param registry The performance registry.
	 */
	void
	initializeReactionNetwork(const options::Options& opts,
		std::shared_ptr<perf::IHandlerRegistry> registry)
	{
		// Get the current process ID
		int procId;
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);

		using NetworkType = core::network::FeReactionNetwork;

		// Get the boundaries from the options
		NetworkType::AmountType maxV = opts.getMaxV();
		NetworkType::AmountType maxI = opts.getMaxI();
		NetworkType::AmountType maxHe = opts.getMaxImpurity();
		NetworkType::AmountType groupingWidthHe = opts.getGroupingWidthA();
		NetworkType::AmountType groupingWidthV = opts.getGroupingWidthB();
		// Take care of the case with no grouping
		if (opts.getGroupingMin() > maxV) {
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

		std::unique_ptr<NetworkType> rNetwork(
			new NetworkType({maxHe, maxV, maxI},
				{{groupingWidthHe, groupingWidthV, maxI + 1}}, 1, opts));
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
	core::network::IReactionNetwork&
	getNetworkHandler() const
	{
		return *theNetworkHandler;
	}
};

} // end namespace network
} // end namespace factory
} // end namespace xolotl

#endif // FEREACTIONHANDLERFACTORY_H
