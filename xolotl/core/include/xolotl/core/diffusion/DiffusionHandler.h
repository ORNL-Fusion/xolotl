#pragma once

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

	//! Device copy of diffusingClusters
	Kokkos::View<IdType*> diffClusterIds;

	//! Device clusters
	using DeviceCluster = network::ClusterCommon<plsm::DeviceMemSpace>;
	Kokkos::View<DeviceCluster*> diffClusters;

	//! Migration energy threshold
	double migrationThreshold;

protected:
	void
	syncDiffusingClusters(network::IReactionNetwork& network);

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
	void
	initialize(network::IReactionNetwork& network,
		std::vector<core::RowColPair>& idPairs) override;

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
