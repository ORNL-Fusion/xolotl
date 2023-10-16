// Includes
#include <xolotl/core/diffusion/Diffusion2DHandler.h>

namespace xolotl
{
namespace core
{
namespace diffusion
{
void
Diffusion2DHandler::syncDiffusionGrid()
{
	diffusGrid = Kokkos::View<int***>("Diffusion Grid", diffusionGrid.size(),
		diffusionGrid[0].size(), diffusingClusters.size());
	auto diffGrid_h = create_mirror_view(diffusGrid);
	for (IdType j = 0; j < diffusionGrid.size(); ++j) {
		for (IdType i = 0; i < diffusionGrid[j].size(); ++i) {
			for (IdType n = 0; n < diffusingClusters.size(); ++n) {
				diffGrid_h(j, i, n) = diffusionGrid[j][i][n];
			}
		}
	}
	deep_copy(diffusGrid, diffGrid_h);
}

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

	syncDiffusionGrid();
}

void
Diffusion2DHandler::computeDiffusion(network::IReactionNetwork& network,
	const StencilConcArray& concVector, Kokkos::View<double*> updatedConcOffset,
	double hxLeft, double hxRight, int ix, double sy, int iy, double sz,
	int) const
{
	// Consider each diffusing cluster.
	// TODO Maintaining a separate index assumes that diffusingClusters is
	// visited in same order as diffusionGrid array for given point.
	// Currently true with C++11, but we'd like to be able to visit the
	// diffusing clusters in any order (so that we can parallelize).
	// Maybe with a zip? or a std::transform?

	if (concVector.size() != 5) {
		throw std::runtime_error(
			"Wrong size for 2D concentration stencil; should be 5, got " +
			std::to_string(concVector.size()));
	}
	Kokkos::Array<Kokkos::View<const double*>, 5> concVec = {concVector[0],
		concVector[1], concVector[2], concVector[3], concVector[4]};

	auto diffGrid = diffusGrid;
	auto clusterIds = this->diffClusterIds;
	auto clusters = this->diffClusters;
	Kokkos::parallel_for(
		clusterIds.size(), KOKKOS_LAMBDA(IdType i) {
			auto currId = clusterIds[i];
			auto cluster = clusters[i];

			// Get the initial concentrations
			double oldConc =
				concVec[0][currId] * diffGrid(iy + 1, ix + 1, i); // middle
			double oldLeftConc =
				concVec[1][currId] * diffGrid(iy + 1, ix, i); // left
			double oldRightConc =
				concVec[2][currId] * diffGrid(iy + 1, ix + 2, i); // right
			double oldBottomConc =
				concVec[3][currId] * diffGrid(iy, ix + 1, i); // bottom
			double oldTopConc =
				concVec[4][currId] * diffGrid(iy + 2, ix + 1, i); // top

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
		});
}

void
Diffusion2DHandler::computePartialsForDiffusion(
	network::IReactionNetwork& network, Kokkos::View<double*> val,
	double hxLeft, double hxRight, int ix, double sy, int iy, double, int) const
{
	auto diffGrid = diffusGrid;
	auto clusterIds = this->diffClusterIds;
	auto clusters = this->diffClusters;

	Kokkos::parallel_for(
		clusterIds.size(), KOKKOS_LAMBDA(IdType i) {
			auto cluster = clusters[i];

			auto leftDiff = cluster.getDiffusionCoefficient(ix);
			auto midDiff = cluster.getDiffusionCoefficient(ix + 1);
			auto rightDiff = cluster.getDiffusionCoefficient(ix + 2);

			// Compute the partial derivatives for diffusion of this cluster
			// for the middle, left, right, bottom, and top grid point
			val[i * 5] = -2.0 * midDiff * ((1.0 / (hxLeft * hxRight)) + sy) *
				diffGrid(iy + 1, ix + 1, i); // middle
			val[(i * 5) + 1] =
				(midDiff * 2.0 / (hxLeft * (hxLeft + hxRight)) +
					(leftDiff - rightDiff) /
						((hxLeft + hxRight) * (hxLeft + hxRight))) *
				diffGrid(iy + 1, ix, i); // left
			val[(i * 5) + 2] =
				(midDiff * 2.0 / (hxRight * (hxLeft + hxRight)) +
					(rightDiff - leftDiff) /
						((hxLeft + hxRight) * (hxLeft + hxRight))) *
				diffGrid(iy + 1, ix + 2, i); // right
			val[(i * 5) + 3] = midDiff * sy * diffGrid(iy, ix + 1, i); // bottom
			val[(i * 5) + 4] =
				midDiff * sy * diffGrid(iy + 2, ix + 1, i); // top
		});
}

} /* end namespace diffusion */
} /* end namespace core */
} /* end namespace xolotl */
