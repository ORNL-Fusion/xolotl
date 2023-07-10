#ifndef SORETDIFFUSIONHANDLER_H
#define SORETDIFFUSIONHANDLER_H

// Includes
#include <xolotl/core/modified/ISoretDiffusionHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
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
	std::vector<double> beta;

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
	virtual void
	initialize(network::IReactionNetwork& network,
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

		using NetworkType = network::IPSIReactionNetwork;
		using AmountType = NetworkType::AmountType;

		auto psiNetwork = dynamic_cast<NetworkType*>(&network);
		auto numSpecies = psiNetwork->getSpeciesListSize();
		auto specIdHe = psiNetwork->getHeliumSpeciesId();

		// Initialize the composition
		auto comp = std::vector<AmountType>(numSpecies, 0);

		// Helium
		comp[specIdHe()] = 1;
		auto clusterId = psiNetwork->findClusterId(comp);
		// Check that the helium cluster is present in the network
		if (clusterId != NetworkType::invalidIndex()) {
			auto cluster = network.getClusterCommon(clusterId);

			// Get its diffusion factor and migration energy
			double diffFactor = cluster.getDiffusionFactor();

			// Don't do anything if the diffusion factor is 0.0
			if (not util::equal(diffFactor, 0.0)) {
				// Note that cluster is diffusing.
				diffusingClusters.emplace_back(clusterId);
				beta.emplace_back(0.0065);

				// This cluster interacts with temperature now
				dfill[clusterId].emplace_back(dof);
				ofill[clusterId].emplace_back(dof);
			}
		}

		// Hydrogen
		comp[specIdHe()] = 0;
		if (psiNetwork->hasDeuterium()) {
			auto clusterSpecies = network.parseSpeciesId("D");
			// Get the cluster
			comp[clusterSpecies()] = 1;
			auto clusterId = psiNetwork->findClusterId(comp);
			// Check that the deuterium cluster is present in the network
			if (clusterId != NetworkType::invalidIndex()) {
				auto cluster = network.getClusterCommon(clusterId);

				// Get its diffusion factor and migration energy
				double diffFactor = cluster.getDiffusionFactor();

				// Don't do anything if the diffusion factor is 0.0
				if (not util::equal(diffFactor, 0.0)) {
					// Note that cluster is diffusing.
					diffusingClusters.emplace_back(clusterId);
					beta.emplace_back(0.0045);

					// This cluster interacts with temperature now
					dfill[clusterId].emplace_back(dof);
					ofill[clusterId].emplace_back(dof);
				}
			}
			comp[clusterSpecies()] = 0;
		}
		if (psiNetwork->hasTritium()) {
			auto clusterSpecies = network.parseSpeciesId("T");
			// Get the cluster
			comp[clusterSpecies()] = 1;
			auto clusterId = psiNetwork->findClusterId(comp);
			// Check that the deuterium cluster is present in the network
			if (clusterId != NetworkType::invalidIndex()) {
				auto cluster = network.getClusterCommon(clusterId);

				// Get its diffusion factor and migration energy
				double diffFactor = cluster.getDiffusionFactor();

				// Don't do anything if the diffusion factor is 0.0
				if (not util::equal(diffFactor, 0.0)) {
					// Note that cluster is diffusing.
					diffusingClusters.emplace_back(clusterId);
					beta.emplace_back(0.0045);

					// This cluster interacts with temperature now
					dfill[clusterId].emplace_back(dof);
					ofill[clusterId].emplace_back(dof);
				}
			}
			comp[clusterSpecies()] = 0;
		}

		// Self-interstitial
		auto specIdI = psiNetwork->getInterstitialSpeciesId();
		comp[specIdI()] = 1;
		clusterId = psiNetwork->findClusterId(comp);
		// Check that the interstitial cluster is present in the network
		if (clusterId != NetworkType::invalidIndex()) {
			auto cluster = network.getClusterCommon(clusterId);

			// Get its diffusion factor and migration energy
			double diffFactor = cluster.getDiffusionFactor();

			// Don't do anything if the diffusion factor is 0.0
			if (not util::equal(diffFactor, 0.0)) {
				// Note that cluster is diffusing.
				diffusingClusters.emplace_back(clusterId);
				beta.emplace_back(0.0128);

				// This cluster interacts with temperature now
				dfill[clusterId].emplace_back(dof);
				ofill[clusterId].emplace_back(dof);
			}
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
