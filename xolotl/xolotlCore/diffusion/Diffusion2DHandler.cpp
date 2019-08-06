// Includes
#include "Diffusion2DHandler.h"

namespace xolotlCore {

void Diffusion2DHandler::initializeDiffusionGrid(
		std::vector<IAdvectionHandler *> advectionHandlers,
		std::vector<double> grid, int ny, double hy, int nz, double hz) {
	// Get the number of diffusing clusters
	int nDiff = diffusingClusters.size();

	// Get the size of the grid in the depth direction
	int nx = grid.size();

	// Initialize the diffusion grid with true everywhere
	diffusionGrid.clear();
	// Initialize it to True
	for (int j = 0; j < ny + 2; j++) {
		std::vector<std::vector<bool> > tempGridBis;
		for (int i = 0; i < nx; i++) {
			tempGridBis.emplace_back(nDiff, true);
		}
		diffusionGrid.push_back(tempGridBis);
	}

	// Initialize the grid position
	NDPoint<3> gridPosition { 0.0, 0.0, 0.0 };

	// Loop on the advection handlers
	for (auto const& currAdvectionHandler : advectionHandlers) {
		// Get the list of advecting clusters
		auto const& advecClusters =
				currAdvectionHandler->getAdvectingClusters();

		// Loop on the spatial grid
		for (int j = -1; j < ny + 1; j++) {
			// Set the grid position
			gridPosition[1] = hy * (double) j;
			for (int i = 0; i < nx; i++) {
				// Set the grid position
				gridPosition[0] = grid[i] - grid[1];

				// Check if we are on a sink
				if (currAdvectionHandler->isPointOnSink(gridPosition)) {
					// We have to find the corresponding reactant in the 
					// diffusion cluster collection.
					for (IReactant const& currAdvCluster : advecClusters) {
						// Initialize n the index in the diffusion index vector
						// TODO use std::find or std::find_if
						int n = 0;
						while (n < nDiff) {
							IReactant const& currDiffCluster =
									diffusingClusters[n];
							if (&currDiffCluster == &currAdvCluster) {
								break;
							}
							n++;
						}
						// Set this diffusion grid value to false
						diffusionGrid[j + 1][i][n] = false;
					}
				}
			}
		}
	}

	return;
}

void Diffusion2DHandler::computeDiffusion(const IReactionNetwork& network,
		double **concVector, double *updatedConcOffset, double hxLeft,
		double hxRight, int ix, int xs, double sy, int iy, double, int) const {

	// Consider each diffusing cluster.
	// TODO Maintaining a separate index assumes that diffusingClusters is
	// visited in same order as diffusionGrid array for given point.
	// Currently true with C++11, but we'd like to be able to visit the
	// diffusing clusters in any order (so that we can parallelize).
	// Maybe with a zip? or a std::transform?
	int diffClusterIdx = 0;
	for (IReactant const& currReactant : diffusingClusters) {

		auto const& cluster = static_cast<IReactant const&>(currReactant);
		int index = cluster.getId() - 1;

		// Get the initial concentrations
		double oldConc = concVector[0][index]
				* diffusionGrid[iy + 1][ix + 1][diffClusterIdx]; // middle
		double oldLeftConc = concVector[1][index]
				* diffusionGrid[iy + 1][ix][diffClusterIdx]; // left
		double oldRightConc = concVector[2][index]
				* diffusionGrid[iy + 1][ix + 2][diffClusterIdx]; // right
		double oldBottomConc = concVector[3][index]
				* diffusionGrid[iy][ix + 1][diffClusterIdx]; // bottom
		double oldTopConc = concVector[4][index]
				* diffusionGrid[iy + 2][ix + 1][diffClusterIdx]; // top

		// Use a simple midpoint stencil to compute the concentration
		double conc = cluster.getDiffusionCoefficient(ix + 1 - xs)
				* (2.0
						* (oldLeftConc + (hxLeft / hxRight) * oldRightConc
								- (1.0 + (hxLeft / hxRight)) * oldConc)
						/ (hxLeft * (hxLeft + hxRight))
						+ sy * (oldBottomConc + oldTopConc - 2.0 * oldConc))
				+ ((cluster.getDiffusionCoefficient(ix + 2 - xs)
						- cluster.getDiffusionCoefficient(ix - xs))
						* (oldRightConc - oldLeftConc)
						/ ((hxLeft + hxRight) * (hxLeft + hxRight)));

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;

		++diffClusterIdx;
	}

	return;
}

void Diffusion2DHandler::computePartialsForDiffusion(
		const IReactionNetwork& network, double *val, int *indices,
		double hxLeft, double hxRight, int ix, int xs, double sy, int iy,
		double, int) const {

	// Consider each diffusing cluster.
	// TODO Maintaining a separate index assumes that diffusingClusters is
	// visited in same order as diffusionGrid array for given point.
	// Currently true with C++11, but we'd like to be able to visit the
	// diffusing clusters in any order (so that we can parallelize).
	// Maybe with a zip? or a std::transform?
	int diffClusterIdx = 0;
	for (IReactant const& currReactant : diffusingClusters) {

		auto const& cluster = static_cast<IReactant const&>(currReactant);
		int index = cluster.getId() - 1;

		// Set the cluster index, the PetscSolver will use it to compute
		// the row and column indices for the Jacobian
		indices[diffClusterIdx] = index;

		// Compute the partial derivatives for diffusion of this cluster
		// for the middle, left, right, bottom, and top grid point
		val[diffClusterIdx * 5] = -2.0
				* cluster.getDiffusionCoefficient(ix + 1 - xs)
				* ((1.0 / (hxLeft * hxRight)) + sy)
				* diffusionGrid[iy + 1][ix + 1][diffClusterIdx]; // middle
		val[(diffClusterIdx * 5) + 1] = (cluster.getDiffusionCoefficient(
				ix + 1 - xs) * 2.0 / (hxLeft * (hxLeft + hxRight))
				+ (cluster.getDiffusionCoefficient(ix - xs)
						- cluster.getDiffusionCoefficient(ix + 2 - xs))
						/ ((hxLeft + hxRight) * (hxLeft + hxRight)))
				* diffusionGrid[iy + 1][ix][diffClusterIdx]; // left
		val[(diffClusterIdx * 5) + 2] = (cluster.getDiffusionCoefficient(
				ix + 1 - xs) * 2.0 / (hxRight * (hxLeft + hxRight))
				+ (cluster.getDiffusionCoefficient(ix + 2 - xs)
						- cluster.getDiffusionCoefficient(ix - xs))
						/ ((hxLeft + hxRight) * (hxLeft + hxRight)))
				* diffusionGrid[iy + 1][ix + 2][diffClusterIdx]; // right
		val[(diffClusterIdx * 5) + 3] = cluster.getDiffusionCoefficient(
				ix + 1 - xs) * sy * diffusionGrid[iy][ix + 1][diffClusterIdx]; // bottom
		val[(diffClusterIdx * 5) + 4] = cluster.getDiffusionCoefficient(
				ix + 1 - xs) * sy
				* diffusionGrid[iy + 2][ix + 1][diffClusterIdx]; // top

		// Increase the index
		diffClusterIdx++;
	}

	return;
}

}/* end namespace xolotlCore */
