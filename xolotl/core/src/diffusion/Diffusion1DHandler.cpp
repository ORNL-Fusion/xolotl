// Includes
#include <iostream>
#include <array>
#include <xolotl/core/diffusion/Diffusion1DHandler.h>

namespace xolotl {
namespace core {
namespace diffusion {

void Diffusion1DHandler::initializeDiffusionGrid(
		std::vector<advection::IAdvectionHandler *> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny, double hy, int ys,
		int nz, double hz, int zs) {
	// Get the number of diffusing clusters
	const int nDiff = diffusingClusters.size();

	// Initialize the diffusion grid with true everywhere
	diffusionGrid.clear();
	for (int i = 0; i < nx + 2; i++) {
		diffusionGrid.emplace_back(nDiff, true);
	}

	// Initialize the grid position
    util::Point<3> gridPosition { 0.0, 0.0, 0.0 };

	// Consider each advection handler
	for (auto const& currAdvectionHandler : advectionHandlers) {

		// Access collection of advecting clusters.
		auto const& advecClusters =
				currAdvectionHandler->getAdvectingClusters();

		// Loop on the spatial grid
		for (int i = 0; i < nx; i++) {
			// Set the grid position
			gridPosition[0] = (grid[i + xs] + grid[i + xs + 1]) / 2.0 - grid[1];

			// Check if we are on a sink
			if (currAdvectionHandler->isPointOnSink(gridPosition)) {
				// We have to find the corresponding reactant in the diffusion
				// cluster collection.
				for (auto const currAdvCluster : advecClusters) {

					// Initialize n the index in the diffusion index vector
					// TODO use std::find or std::find_if?
					int n = 0;
					while (n < nDiff) {
						auto const currDiffCluster = diffusingClusters[n];
						if (currDiffCluster == currAdvCluster) {
							break;
						}
						n++;
					}
					// Set this diffusion grid value to false
					diffusionGrid[i][n] = false;
				}
			}
		}
	}

	return;
}

void Diffusion1DHandler::computeDiffusion(
		network::IReactionNetwork& network, double **concVector,
		double *updatedConcOffset, double hxLeft, double hxRight, int ix,
		double, int, double, int) const {

	// Consider each diffusing cluster.
	// TODO Maintaining a separate index assumes that diffusingClusters is
	// visited in same order as diffusionGrid array for given point.
	// Currently true with C++11, but we'd like to be able to visit the
	// diffusing clusters in any order (so that we can parallelize).
	// Maybe with a zip? or a std::transform?
	int diffClusterIdx = 0;
	for (auto const& currId : diffusingClusters) {

		auto cluster = network.getClusterCommon(currId);

		// Get the initial concentrations
		double oldConc = concVector[0][currId]
				* diffusionGrid[ix + 1][diffClusterIdx];
		double oldLeftConc = concVector[1][currId]
				* diffusionGrid[ix][diffClusterIdx];
		double oldRightConc = concVector[2][currId]
				* diffusionGrid[ix + 2][diffClusterIdx];

		// Use a simple midpoint stencil to compute the concentration
		// TODO: Check we are using the correct grid index
		double conc = (cluster.getDiffusionCoefficient(ix + 1) * 2.0
				* (oldLeftConc + (hxLeft / hxRight) * oldRightConc
						- (1.0 + (hxLeft / hxRight)) * oldConc)
				/ (hxLeft * (hxLeft + hxRight)))
				+ ((cluster.getDiffusionCoefficient(ix + 2)
						- cluster.getDiffusionCoefficient(ix))
						* (oldRightConc - oldLeftConc)
						/ ((hxLeft + hxRight) * (hxLeft + hxRight)));

		// Update the concentration of the cluster
		updatedConcOffset[currId] += conc;

		// Increase the index
		diffClusterIdx++;
	}

	return;
}

void Diffusion1DHandler::computePartialsForDiffusion(
		network::IReactionNetwork& network, double *val,
		int *indices, double hxLeft, double hxRight, int ix, double, int,
		double, int) const {

	// Loop on them
	// TODO Maintaining a separate index assumes that diffusingClusters is
	// visited in same order as diffusionGrid array for given point.
	// Currently true with C++11, but we'd like to be able to visit the
	// diffusing clusters in any order (so that we can parallelize).
	// Maybe with a zip? or a std::transform?
	int diffClusterIdx = 0;
	for (auto const& currId : diffusingClusters) {

		auto cluster = network.getClusterCommon(currId);

		// Set the cluster index, the PetscSolver will use it to compute
		// the row and column indices for the Jacobian
		indices[diffClusterIdx] = currId;

		// Compute the partial derivatives for diffusion of this cluster
		// for the middle, left, and right grid point
		val[diffClusterIdx * 3] = -2.0 * cluster.getDiffusionCoefficient(ix + 1)
				/ (hxLeft * hxRight) * diffusionGrid[ix + 1][diffClusterIdx]; // middle
		val[(diffClusterIdx * 3) + 1] = (cluster.getDiffusionCoefficient(ix)
				* 2.0 / (hxLeft * (hxLeft + hxRight))
				+ (cluster.getDiffusionCoefficient(ix)
						- cluster.getDiffusionCoefficient(ix + 2))
						/ ((hxLeft + hxRight) * (hxLeft + hxRight)))
				* diffusionGrid[ix][diffClusterIdx]; // left
		val[(diffClusterIdx * 3) + 2] = (cluster.getDiffusionCoefficient(ix)
				* 2.0 / (hxRight * (hxLeft + hxRight))
				+ (cluster.getDiffusionCoefficient(ix + 2)
						- cluster.getDiffusionCoefficient(ix))
						/ ((hxLeft + hxRight) * (hxLeft + hxRight)))
				* diffusionGrid[ix + 2][diffClusterIdx]; // right

		// Increase the index
		diffClusterIdx++;
	}

	return;
}

}/* end namespace diffusion */
}/* end namespace core */
}/* end namespace xolotl */
