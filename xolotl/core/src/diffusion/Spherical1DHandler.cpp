
#include <array>
#include <iostream>

#include <xolotl/core/diffusion/Spherical1DHandler.h>

namespace xolotl
{
namespace core
{
namespace diffusion
{
void
Spherical1DHandler::syncGrids()
{
	// Diffusion grid
	diffusGrid = Kokkos::View<int**>(
		"Diffusion Grid", diffusionGrid.size(), diffusingClusters.size());
	auto diffGrid_h = create_mirror_view(diffusGrid);
	for (IdType i = 0; i < diffusionGrid.size(); ++i) {
		for (IdType n = 0; n < diffusingClusters.size(); ++n) {
			diffGrid_h(i, n) = diffusionGrid[i][n];
		}
	}
	deep_copy(diffusGrid, diffGrid_h);

	// Physical grid
	physGrid = Kokkos::View<double*>("Physical Grid", physicalGrid.size());
	auto physGrid_h = create_mirror_view(physGrid);
	for (IdType i = 0; i < physicalGrid.size(); ++i) {
		physGrid_h(i) = physicalGrid[i];
	}

	deep_copy(physGrid, physGrid_h);
}

void
Spherical1DHandler::initializeDiffusionGrid(
	std::vector<advection::IAdvectionHandler*> advectionHandlers,
	std::vector<double> grid, int nx, int xs, int ny, double hy, int ys, int nz,
	double hz, int zs)
{
	// Physical grid first
	for (int i = 0; i < nx + 2; i++) {
		physicalGrid.push_back(grid[i + xs]);
	}

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

	syncGrids();
}

void
Spherical1DHandler::computeDiffusion(network::IReactionNetwork& network,
	const StencilConcArray& concVector, Kokkos::View<double*> updatedConcOffset,
	double hxLeft, double hxRight, int ix, bool isBC, double sy, int iy,
	double sz, int) const
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

	// Compute alpha
	auto alpha = network.getTritiumFlux(concVector[0][1], concVector[0][2],
					 concVector[0][0], clusters[0].getTemperature(ix + 1)) /
		clusters[0].getDiffusionCoefficient(ix + 1);

	Kokkos::parallel_for(
		clusterIds.size(), KOKKOS_LAMBDA(IdType i) {
			auto currId = clusterIds[i];
			auto cluster = clusters[i];

			// Get the initial concentrations
			double oldConc = concVec[0][currId] * diffGrid(ix + 1, i);
			double oldLeftConc = concVec[1][currId] * diffGrid(ix, i);
			double oldRightConc = concVec[2][currId] * diffGrid(ix + 2, i);
			double midDiff = cluster.getDiffusionCoefficient(ix + 1);
			double leftTemp = cluster.getTemperature(ix);
			double midTemp = cluster.getTemperature(ix + 1);
			double rightTemp = cluster.getTemperature(ix + 2);

			// Surface location
			if (isBC) {
				// Boundary condition depending on the flux at the surface
				double conc =
					(midDiff * 2.0 * (oldLeftConc + hxLeft * alpha - oldConc) /
						(hxLeft * (hxLeft + hxRight)));

				// The one specific to spherical coordinates
				conc += 2.0 * midDiff * (oldConc - oldLeftConc) /
					(physGrid(ix + 1) * hxLeft);

				// Update the concentration of the cluster
				updatedConcOffset[currId] += conc;
			}

			// Everywhere else
			else {
				// Use a simple midpoint stencil to compute the concentration
				// The usual term is there assuming the temperature is the same
				// everywhere
				double conc = (midDiff * 2.0 *
					(oldLeftConc + (hxLeft / hxRight) * oldRightConc -
						(1.0 + (hxLeft / hxRight)) * oldConc) /
					(hxLeft * (hxLeft + hxRight)));

				// The one specific to spherical coordinates
				conc += 2.0 * midDiff * (oldConc - oldLeftConc) /
					(physGrid(ix + 1) * hxLeft);

				// Update the concentration of the cluster
				updatedConcOffset[currId] += conc;
			}
		});
}

void
Spherical1DHandler::computePartialsForDiffusion(
	network::IReactionNetwork& network, Kokkos::View<double*> val,
	double hxLeft, double hxRight, int ix, bool isBC, double, int, double,
	int) const
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

			// Surface location
			if (isBC) {
				val[i * 3] = (-2.0 * midDiff / (hxLeft * (hxLeft + hxRight)) +
								 2.0 * midDiff / (physGrid(ix + 1) * hxLeft)) *
					diffGrid(ix + 1, i); // middle
				val[(i * 3) + 1] =
					(midDiff * 2.0 / (hxLeft * (hxLeft + hxRight)) -
						2.0 * midDiff / (physGrid(ix + 1) * hxLeft)) *
					diffGrid(ix, i); // left
				val[(i * 3) + 2] = 0.0; // right
			}

			// Everywhere else
			else {
				// Compute the partial derivatives for diffusion of this cluster
				// for the middle, left, and right grid point
				val[i * 3] = (-2.0 * midDiff / (hxLeft * hxRight) +
								 2.0 * midDiff / (physGrid(ix + 1) * hxLeft)) *
					diffGrid(ix + 1, i); // middle
				val[(i * 3) + 1] =
					(midDiff * 2.0 / (hxLeft * (hxLeft + hxRight)) -
						2.0 * midDiff / (physGrid(ix + 1) * hxLeft)) *
					diffGrid(ix, i); // left
				val[(i * 3) + 2] =
					(midDiff * 2.0 / (hxRight * (hxLeft + hxRight))) *
					diffGrid(ix + 2, i); // right
			}
		});
}

} /* end namespace diffusion */
} /* end namespace core */
} /* end namespace xolotl */
