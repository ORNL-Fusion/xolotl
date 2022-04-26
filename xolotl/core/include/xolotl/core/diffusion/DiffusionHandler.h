#ifndef DIFFUSIONHANDLER_H
#define DIFFUSIONHANDLER_H

// Includes
#include <xolotl/core/diffusion/IDiffusionHandler.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace diffusion
{
/**
 * This class realizes the IDiffusionHandler interface responsible for all
 * the physical parts for the diffusion of mobile clusters. It needs to have
 * subclasses implementing the compute diffusion methods.
 */
class DiffusionHandler : public IDiffusionHandler
{
protected:
	//! Collection of diffusing clusters.
	std::vector<IdType> diffusingClusters;

	//! Migration energy threshold
	double migrationThreshold;

public:
	//! The Constructor
	DiffusionHandler(double threshold) : migrationThreshold(threshold)
	{
	}

	//! The Destructor
	~DiffusionHandler()
	{
	}

	/**
	 * Initialize the off-diagonal part of the Jacobian. If this step is skipped
	 * it won't be possible to set the partial derivatives for the diffusion.
	 *
	 * The value 1 is set in ofillMap if a cluster has a non zero diffusion
	 * coefficient and a migration energy lower than the threshold.
	 *
	 * \see IDiffusionHandler.h
	 */
	virtual void
	initializeOFill(network::IReactionNetwork& network,
		network::IReactionNetwork::SparseFillMap& ofillMap) override
	{
		// Clear the index vector
		diffusingClusters.clear();

		// Consider each cluster
		for (std::size_t i = 0; i < network.getNumClusters(); i++) {
			auto cluster = network.getClusterCommon(i);

			// Get its diffusion factor and migration energy
			double diffFactor = cluster.getDiffusionFactor();
			double migration = cluster.getMigrationEnergy();

			// Don't do anything if the diffusion factor is 0.0
			if (util::equal(diffFactor, 0.0) || migration > migrationThreshold)
				continue;

			// Note that cluster is diffusing.
			diffusingClusters.emplace_back(i);

			// Set the ofill value to 1 for this cluster
			ofillMap[i].emplace_back(i);
		}

		return;
	}

	/**
	 * Get the total number of diffusing clusters in the network.
	 *
	 * \see IDiffusionHandler.h
	 */
	int
	getNumberOfDiffusing() const override
	{
		return diffusingClusters.size();
	}

	/**
	 * Get the vector of IDs of diffusing clusters in the network.
	 *
	 * \see IDiffusionHandler.h
	 */
	virtual std::vector<IdType>
	getDiffusingIds() const override
	{
		return diffusingClusters;
	}
};
// end class DiffusionHandler

} /* end namespace diffusion */
} /* end namespace core */
} /* end namespace xolotl */
#endif
