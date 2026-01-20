#pragma once

// Includes
#include <xolotl/core/modified/ISoretDiffusionHandler.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace modified
{
/**
 * This class realizes the IDiffusionHandler interface responsible for all
 * the physical parts for the diffusion of mobile clusters. It needs to have
 * subclasses implementing the compute diffusion methods.
 */
class SoretDiffusionHandler : public ISoretDiffusionHandler
{
protected:
	//! Collection of diffusing clusters.
	std::vector<IdType> diffusingClusters;

	//! Device copy of diffusingClusters
	Kokkos::View<IdType*> diffClusterIds;

	//! Device clusters
	using DeviceCluster = network::ClusterCommon<plsm::DeviceMemSpace>;
	Kokkos::View<DeviceCluster*> diffClusters;

	//! The number of degrees of freedom in the network
	IdType dof;

	//! The x grid
	std::vector<double> xGrid;

	//! The surface position
	int surfacePosition;

	//! The local xs
	int localXs;

	//! The beta factor
	std::vector<double> betaFactor;

	//! Device copy of betaFactor
	Kokkos::View<double*> beta;

protected:
	void
	syncDiffusingClusters(network::IReactionNetwork& network);

public:
	//! The Constructor
	SoretDiffusionHandler() : dof(0), surfacePosition(0), localXs(0)
	{
	}

	//! The Destructor
	~SoretDiffusionHandler()
	{
	}

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	void
	initialize(network::IReactionNetwork& network,
		std::vector<core::RowColPair>& idPairs, std::vector<double> grid,
		int xs) override;

	int
	getNumberOfDiffusing() const override
	{
		return diffusingClusters.size();
	}

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	void
	updateSurfacePosition(int surfacePos) override
	{
		surfacePosition = surfacePos;
	}

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	void
	computeDiffusion(network::IReactionNetwork& network,
		const StencilConcArray& concVector,
		Kokkos::View<double*> updatedConcOffset, double hxLeft, double hxRight,
		int ix, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) const override;

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	bool
	computePartialsForDiffusion(network::IReactionNetwork& network,
		const StencilConcArray& concVector, Kokkos::View<double*> values,
		double hxLeft, double hxRight, int ix, double sy = 0.0, int iy = 0,
		double sz = 0.0, int iz = 0) const override;
};
// end class SoretDiffusionHandler

} /* end namespace modified */
} /* end namespace core */
} /* end namespace xolotl */
