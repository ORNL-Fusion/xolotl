// Includes
#include <array>
#include <iostream>

#include <xolotl/core/diffusion/Diffusion1DHandler.h>

namespace xolotl
{
namespace core
{
namespace diffusion
{
void
Diffusion1DHandler::initializeDiffusionGrid(
	std::vector<advection::IAdvectionHandler*> advectionHandlers,
	std::vector<double> grid, int nx, int xs, int ny, double hy, int ys, int nz,
	double hz, int zs)
{
	// Get the number of diffusing clusters
	const int nDiff = diffusingClusters.size();

	// Initialize the diffusion grid with true everywhere
	diffusionGrid.clear();
	for (int i = 0; i < nx + 2; i++) {
		diffusionGrid.emplace_back(nDiff, true);
	}

	// Initialize the grid position
	plsm::SpaceVector<double, 3> gridPosition{0.0, 0.0, 0.0};

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
					auto it = find(diffusingClusters.begin(),
						diffusingClusters.end(), currAdvCluster);
					if (it != diffusingClusters.end()) {
						// Set this diffusion grid value to false
						diffusionGrid[i][(*it)] = false;
					}
					else {
						throw std::runtime_error(
							"\nThe advecting cluster of id: " +
							std::to_string(currAdvCluster) +
							" was not found in the diffusing clusters, cannot "
							"use the diffusion!");
					}
				}
			}
		}
	}

	return;
}

void
Diffusion1DHandler::computeDiffusion(network::IReactionNetwork& network,
	double** concVector, double* updatedConcOffset, double hxLeft,
	double hxRight, int ix, double, int, double, int) const
{
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
		double oldConc =
			concVector[0][currId] * diffusionGrid[ix + 1][diffClusterIdx];
		double oldLeftConc =
			concVector[1][currId] * diffusionGrid[ix][diffClusterIdx];
		double oldRightConc =
			concVector[2][currId] * diffusionGrid[ix + 2][diffClusterIdx];
		double leftDiff = cluster.getDiffusionCoefficient(ix),
			   midDiff = cluster.getDiffusionCoefficient(ix + 1),
			   rightDiff = cluster.getDiffusionCoefficient(ix + 2);
		double leftTemp = cluster.getTemperature(ix),
			   midTemp = cluster.getTemperature(ix + 1),
			   rightTemp = cluster.getTemperature(ix + 2);

		// Use a simple midpoint stencil to compute the concentration
		double conc = (midDiff * 2.0 *
						  (oldLeftConc + (hxLeft / hxRight) * oldRightConc -
							  (1.0 + (hxLeft / hxRight)) * oldConc) /
						  (hxLeft * (hxLeft + hxRight))) +
			((rightDiff - leftDiff) * (oldRightConc - oldLeftConc) /
				((hxLeft + hxRight) * (hxLeft + hxRight)));

		// Update the concentration of the cluster
		updatedConcOffset[currId] += conc;

		// Increase the index
		diffClusterIdx++;
	}

	return;
}

void
Diffusion1DHandler::computePartialsForDiffusion(
	network::IReactionNetwork& network, double* val, IdType* indices,
	double hxLeft, double hxRight, int ix, double, int, double, int) const
{
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
		double leftDiff = cluster.getDiffusionCoefficient(ix),
			   midDiff = cluster.getDiffusionCoefficient(ix + 1),
			   rightDiff = cluster.getDiffusionCoefficient(ix + 2);
		double leftTemp = cluster.getTemperature(ix),
			   midTemp = cluster.getTemperature(ix + 1),
			   rightTemp = cluster.getTemperature(ix + 2);

		// Compute the partial derivatives for diffusion of this cluster
		// for the middle, left, and right grid point
		val[diffClusterIdx * 3] = (-2.0 * midDiff / (hxLeft * hxRight)) *
			diffusionGrid[ix + 1][diffClusterIdx]; // middle
		val[(diffClusterIdx * 3) + 1] =
			(midDiff * 2.0 / (hxLeft * (hxLeft + hxRight)) +
				(leftDiff - rightDiff) /
					((hxLeft + hxRight) * (hxLeft + hxRight))) *
			diffusionGrid[ix][diffClusterIdx]; // left
		val[(diffClusterIdx * 3) + 2] =
			(midDiff * 2.0 / (hxRight * (hxLeft + hxRight)) +
				(rightDiff - leftDiff) /
					((hxLeft + hxRight) * (hxLeft + hxRight))) *
			diffusionGrid[ix + 2][diffClusterIdx]; // right

		// Increase the index
		diffClusterIdx++;
	}

	return;
}

} /* end namespace diffusion */
} /* end namespace core */
} /* end namespace xolotl */
