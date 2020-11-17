// Includes
#include <xolotl/core/diffusion/Diffusion2DHandler.h>

namespace xolotl
{
namespace core
{
namespace diffusion
{
void
Diffusion2DHandler::initializeDiffusionGrid(
	std::vector<advection::IAdvectionHandler*> advectionHandlers,
	std::vector<double> grid, int nx, int xs, int ny, double hy, int ys, int nz,
	double hz, int zs)
{
	// Get the number of diffusing clusters
	int nDiff = diffusingClusters.size();

	// Initialize the diffusion grid with true everywhere
	diffusionGrid.clear();
	// Initialize it to True
	for (int j = 0; j < ny + 2; j++) {
		std::vector<std::vector<bool>> tempGridBis;
		for (int i = 0; i < nx + 2; i++) {
			tempGridBis.emplace_back(nDiff, true);
		}
		diffusionGrid.push_back(tempGridBis);
	}

	// Initialize the grid position
	plsm::SpaceVector<double, 3> gridPosition{0.0, 0.0, 0.0};

	// Loop on the advection handlers
	for (auto const& currAdvectionHandler : advectionHandlers) {
		// Get the list of advecting clusters
		auto const& advecClusters =
			currAdvectionHandler->getAdvectingClusters();

		// Loop on the spatial grid
		for (int j = -1; j < ny + 1; j++) {
			// Set the grid position
			gridPosition[1] = hy * (double)(j + ys);
			for (int i = 0; i < nx; i++) {
				// Set the grid position
				gridPosition[0] =
					(grid[i + xs] + grid[i + xs + 1]) / 2.0 - grid[1];

				// Check if we are on a sink
				if (currAdvectionHandler->isPointOnSink(gridPosition)) {
					// We have to find the corresponding reactant in the
					// diffusion cluster collection.
					for (auto const& currAdvCluster : advecClusters) {
						auto it = find(diffusingClusters.begin(),
							diffusingClusters.end(), currAdvCluster);
						if (it != diffusingClusters.end()) {
							// Set this diffusion grid value to false
							diffusionGrid[j + 1][i][(*it)] = false;
						}
						else {
							throw std::runtime_error(
								"\nThe advecting cluster of id: " +
								std::to_string(currAdvCluster) +
								" was not found in the diffusing clusters, "
								"cannot use the diffusion!");
						}
					}
				}
			}
		}
	}

	return;
}

void
Diffusion2DHandler::computeDiffusion(network::IReactionNetwork& network,
	double** concVector, double* updatedConcOffset, double hxLeft,
	double hxRight, int ix, double sy, int iy, double, int) const
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
		double oldConc = concVector[0][currId] *
			diffusionGrid[iy + 1][ix + 1][diffClusterIdx]; // middle
		double oldLeftConc = concVector[1][currId] *
			diffusionGrid[iy + 1][ix][diffClusterIdx]; // left
		double oldRightConc = concVector[2][currId] *
			diffusionGrid[iy + 1][ix + 2][diffClusterIdx]; // right
		double oldBottomConc = concVector[3][currId] *
			diffusionGrid[iy][ix + 1][diffClusterIdx]; // bottom
		double oldTopConc = concVector[4][currId] *
			diffusionGrid[iy + 2][ix + 1][diffClusterIdx]; // top

		// Use a simple midpoint stencil to compute the concentration
		double conc = cluster.getDiffusionCoefficient(ix + 1) *
				(2.0 *
						(oldLeftConc + (hxLeft / hxRight) * oldRightConc -
							(1.0 + (hxLeft / hxRight)) * oldConc) /
						(hxLeft * (hxLeft + hxRight)) +
					sy * (oldBottomConc + oldTopConc - 2.0 * oldConc)) +
			((cluster.getDiffusionCoefficient(ix + 2) -
				 cluster.getDiffusionCoefficient(ix)) *
				(oldRightConc - oldLeftConc) /
				((hxLeft + hxRight) * (hxLeft + hxRight)));

		// Update the concentration of the cluster
		updatedConcOffset[currId] += conc;

		++diffClusterIdx;
	}

	return;
}

void
Diffusion2DHandler::computePartialsForDiffusion(
	network::IReactionNetwork& network, double* val, int* indices,
	double hxLeft, double hxRight, int ix, double sy, int iy, double, int) const
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

		// Set the cluster index, the PetscSolver will use it to compute
		// the row and column indices for the Jacobian
		indices[diffClusterIdx] = currId;

		// Compute the partial derivatives for diffusion of this cluster
		// for the middle, left, right, bottom, and top grid point
		val[diffClusterIdx * 5] = -2.0 *
			cluster.getDiffusionCoefficient(ix + 1) *
			((1.0 / (hxLeft * hxRight)) + sy) *
			diffusionGrid[iy + 1][ix + 1][diffClusterIdx]; // middle
		val[(diffClusterIdx * 5) + 1] =
			(cluster.getDiffusionCoefficient(ix + 1) * 2.0 /
					(hxLeft * (hxLeft + hxRight)) +
				(cluster.getDiffusionCoefficient(ix) -
					cluster.getDiffusionCoefficient(ix + 2)) /
					((hxLeft + hxRight) * (hxLeft + hxRight))) *
			diffusionGrid[iy + 1][ix][diffClusterIdx]; // left
		val[(diffClusterIdx * 5) + 2] =
			(cluster.getDiffusionCoefficient(ix + 1) * 2.0 /
					(hxRight * (hxLeft + hxRight)) +
				(cluster.getDiffusionCoefficient(ix + 2) -
					cluster.getDiffusionCoefficient(ix)) /
					((hxLeft + hxRight) * (hxLeft + hxRight))) *
			diffusionGrid[iy + 1][ix + 2][diffClusterIdx]; // right
		val[(diffClusterIdx * 5) + 3] =
			cluster.getDiffusionCoefficient(ix + 1) * sy *
			diffusionGrid[iy][ix + 1][diffClusterIdx]; // bottom
		val[(diffClusterIdx * 5) + 4] =
			cluster.getDiffusionCoefficient(ix + 1) * sy *
			diffusionGrid[iy + 2][ix + 1][diffClusterIdx]; // top

		// Increase the index
		diffClusterIdx++;
	}

	return;
}

} /* end namespace diffusion */
} /* end namespace core */
} /* end namespace xolotl */
