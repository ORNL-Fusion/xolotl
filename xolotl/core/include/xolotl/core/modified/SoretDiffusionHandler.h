#ifndef SORETDIFFUSIONHANDLER_H
#define SORETDIFFUSIONHANDLER_H

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
	std::vector<std::size_t> diffusingClusters;

	/**
	 * The number of degrees of freedom in the network
	 */
	int dof;

	/**
	 * The x grid
	 */
	std::vector<double> xGrid;

	/**
	 * The surface position
	 */
	int surfacePosition;

	/**
	 * The local xs
	 */
	int localXs;

	/**
	 * The beta factor
	 */
	double beta;

public:
	//! The Constructor
	SoretDiffusionHandler() :
		dof(0),
		surfacePosition(0),
		beta(0.0065),
		//		beta(0.000065),
		//		beta(0.0),
		localXs(0)
	{
	}

	//! The Destructor
	~SoretDiffusionHandler()
	{
	}

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	virtual void
	initialize(const network::IReactionNetwork& network,
		network::IReactionNetwork::SparseFillMap& ofill,
		network::IReactionNetwork::SparseFillMap& dfill,
		std::vector<double> grid, int xs) override
	{
		// Clear the index vector
		diffusingClusters.clear();

		// Set dof
		dof = network.getDOF();

		// Copy the grid
		xGrid = grid;
		localXs = xs;

		// Consider each cluster
		for (std::size_t i = 0; i < network.getNumClusters(); i++) {
			auto cluster = network.getClusterCommon(i);

			// Get its diffusion factor and migration energy
			double diffFactor = cluster.getDiffusionFactor();

			// Don't do anything if the diffusion factor is 0.0
			if (util::equal(diffFactor, 0.0))
				continue;

			// Note that cluster is diffusing.
			diffusingClusters.emplace_back(i);

			// This cluster interacts with temperature now
			dfill[i].emplace_back(dof);
			ofill[i].emplace_back(dof);
		}

		return;
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
	virtual void
	computeDiffusion(network::IReactionNetwork& network, double** concVector,
		double* updatedConcOffset, double hxLeft, double hxRight, int ix,
		double sy = 0.0, int iy = 0, double sz = 0.0, int iz = 0) const;

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	virtual bool
	computePartialsForDiffusion(network::IReactionNetwork& network,
		double** concVector, double* val, IdType* indices, double hxLeft,
		double hxRight, int ix, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) const;
};
// end class SoretDiffusionHandler

} /* end namespace modified */
} /* end namespace core */
} /* end namespace xolotl */
#endif
