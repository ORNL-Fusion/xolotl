
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
Diffusion1DHandler::syncDiffusionGrid()
{
	diffusGrid = Kokkos::View<int**>(
		"Diffusion Grid", diffusionGrid.size(), diffusingClusters.size());
	auto diffGrid_h = create_mirror_view(diffusGrid);
	for (IdType i = 0; i < diffusionGrid.size(); ++i) {
		for (IdType n = 0; n < diffusingClusters.size(); ++n) {
			diffGrid_h(i, n) = diffusionGrid[i][n];
		}
	}
	deep_copy(diffusGrid, diffGrid_h);
}

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

	syncDiffusionGrid();
}

void
Diffusion1DHandler::computeDiffusion(network::IReactionNetwork& network,
	const StencilConcArray& concVector, Kokkos::View<double*> updatedConcOffset,
	double hxLeft, double hxRight, int ix, double sy, int iy, double sz,
	int) const
{
	if (concVector.size() != 3) {
		throw std::runtime_error(
			"Wrong size for 1D concentration stencil; should be 3, got " +
			std::to_string(concVector.size()));
	}
	Kokkos::Array<Kokkos::View<const double*>, 3> concVec = {
		concVector[0], concVector[1], concVector[2]};

	auto diffGrid = diffusGrid;
	auto clusterIds = this->diffClusterIds;
	auto clusters = this->diffClusters;
	Kokkos::parallel_for(
		clusterIds.size(), KOKKOS_LAMBDA(IdType i) {
			auto currId = clusterIds[i];
			auto cluster = clusters[i];

			// Get the initial concentrations
			double oldConc = concVec[0][currId] * diffGrid(ix + 1, i);
			double oldLeftConc = concVec[1][currId] * diffGrid(ix, i);
			double oldRightConc = concVec[2][currId] * diffGrid(ix + 2, i);
			double leftDiff = cluster.getDiffusionCoefficient(ix);
			double midDiff = cluster.getDiffusionCoefficient(ix + 1);
			double rightDiff = cluster.getDiffusionCoefficient(ix + 2);
			double leftTemp = cluster.getTemperature(ix);
			double midTemp = cluster.getTemperature(ix + 1);
			double rightTemp = cluster.getTemperature(ix + 2);

			// Use a simple midpoint stencil to compute the concentration
			double conc = (midDiff * 2.0 *
							  (oldLeftConc + (hxLeft / hxRight) * oldRightConc -
								  (1.0 + (hxLeft / hxRight)) * oldConc) /
							  (hxLeft * (hxLeft + hxRight))) +
				((rightDiff - leftDiff) * (oldRightConc - oldLeftConc) /
					((hxLeft + hxRight) * (hxLeft + hxRight)));

			// Update the concentration of the cluster
			// TODO: Should this use atomic_?
			updatedConcOffset[currId] += conc;
		});

	// TODO: Maybe we need a Kokkos::fence() here?
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
