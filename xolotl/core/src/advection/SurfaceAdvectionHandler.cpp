// Includes
#include <xolotl/core/advection/SurfaceAdvectionHandler.h>

namespace xolotl
{
namespace core
{
namespace advection
{
void
SurfaceAdvectionHandler::initializeAdvectionGrid(
	std::vector<IAdvectionHandler*> advectionHandlers, std::vector<double> grid,
	int nx, int xs, int ny, double hy, int ys, int nz, double hz, int zs)
{
	// Get the number of advecting clusters
	int nAdvec = advectingClusters.size();

	// Initialize the advection grid with true everywhere
	advectionGrid.clear();
	// Initialize it to True
	for (int k = 0; k < nz + 2; k++) {
		std::vector<std::vector<std::vector<bool>>> tempGridTer;
		for (int j = 0; j < ny + 2; j++) {
			std::vector<std::vector<bool>> tempGridBis;
			for (int i = 0; i < nx + 2; i++) {
				tempGridBis.emplace_back(nAdvec, true);
			}
			tempGridTer.push_back(tempGridBis);
		}
		advectionGrid.push_back(tempGridTer);
	}

	// Initialize the grid position
	plsm::SpaceVector<double, 3> gridPosition{0.0, 0.0, 0.0};

	// Consider each advection handler.
	for (auto const& currAdvecHandler : advectionHandlers) {
		// Get the list of advecting clusters
		auto const& otherAdvecClusters =
			currAdvecHandler->getAdvectingClusters();

		// Loop on the spatial grid
		for (int k = -1; k < nz + 1; k++) {
			// Set the grid position
			gridPosition[2] = hz * (double)(k + zs);
			for (int j = -1; j < ny + 1; j++) {
				// Set the grid position
				gridPosition[1] = hy * (double)(j + ys);
				for (int i = 0; i < nx + 2; i++) {
					// Set the grid position
					if (i + xs == nx - 1)
						gridPosition[0] = grid[i + xs] +
							(grid[i + xs] - grid[i + xs - 1]) / 2.0;
					else
						gridPosition[0] =
							(grid[i + xs] + grid[i + xs + 1]) / 2.0;

					// Check if we are on a sink
					if (currAdvecHandler->isPointOnSink(gridPosition)) {
						// We have to find the corresponding index in the
						// diffusion index vector
						for (auto const currAdvCluster : otherAdvecClusters) {
							auto it = find(advectingClusters.begin(),
								advectingClusters.end(), currAdvCluster);
							if (it != advectingClusters.end()) {
								// Set this diffusion grid value to false
								advectionGrid[k + 1][j + 1][i + 1][(*it)] =
									false;
							}
							else {
								throw std::runtime_error(
									"\nThe advecting cluster of id: " +
									std::to_string(currAdvCluster) +
									" was not found in the advecting clusters, "
									"cannot "
									"use the advection!");
							}
						}
					}
				}
			}
		}
	}

	return;
}

void
SurfaceAdvectionHandler::computeAdvection(network::IReactionNetwork& network,
	const plsm::SpaceVector<double, 3>& pos, double** concVector,
	double* updatedConcOffset, double hxLeft, double hxRight, int ix, double hy,
	int iy, double hz, int iz) const
{
	// Consider each advecting cluster
	// TODO Maintaining a separate index assumes that advectingClusters is
	// visited in same order as advectionGrid array for given point
	// and the sinkStrengthVector.
	// Currently true with C++11, but we'd like to be able to visit the
	// advecting clusters in any order (so that we can parallelize).
	// Maybe with a zip? or a std::transform?
	int advClusterIdx = 0;
	for (auto const& currId : advectingClusters) {
		auto cluster = network.getClusterCommon(currId);

		// Get the initial concentrations
		double oldConc = concVector[0][currId] *
			advectionGrid[iz + 1][iy + 1][ix + 1][advClusterIdx]; // middle
		double oldRightConc = concVector[2][currId] *
			advectionGrid[iz + 1][iy + 1][ix + 2][advClusterIdx]; // right

		// Compute the concentration as explained in the description of the
		// method
		double conc = (3.0 * sinkStrengthVector[advClusterIdx] *
						  cluster.getDiffusionCoefficient(ix + 1)) *
			((oldRightConc / pow(pos[0] - location + hxRight, 4)) -
				(oldConc / pow(pos[0] - location, 4))) /
			(kBoltzmann * cluster.getTemperature(ix + 1) * hxRight);

		conc += (3.0 * sinkStrengthVector[advClusterIdx] * oldConc) *
			(cluster.getDiffusionCoefficient(ix + 2) /
					cluster.getTemperature(ix + 2) -
				cluster.getDiffusionCoefficient(ix + 1) /
					cluster.getTemperature(ix + 1)) /
			(kBoltzmann * hxRight * pow(pos[0] - location, 4));

		// Update the concentration of the cluster
		updatedConcOffset[currId] += conc;

		++advClusterIdx;
	}

	return;
}

void
SurfaceAdvectionHandler::computePartialsForAdvection(
	network::IReactionNetwork& network, double* val, IdType* indices,
	const plsm::SpaceVector<double, 3>& pos, double hxLeft, double hxRight,
	int ix, double hy, int iy, double hz, int iz) const
{
	// Get the number of advecting cluster
	int nAdvec = advectingClusters.size();

	// Consider each advecting cluster.
	// TODO Maintaining a separate index assumes that advectingClusters is
	// visited in same order as advectionGrid array for given point
	// and the sinkStrengthVector.
	// Currently true with C++11, but we'd like to be able to visit the
	// advecting clusters in any order (so that we can parallelize).
	// Maybe with a zip? or a std::transform?
	int advClusterIdx = 0;
	for (auto const& currId : advectingClusters) {
		auto cluster = network.getClusterCommon(currId);
		// Get the diffusion coefficient of the cluster
		double diffCoeff = cluster.getDiffusionCoefficient(ix + 1);
		// Get the sink strength value
		double sinkStrength = sinkStrengthVector[advClusterIdx];

		// Set the cluster index that will be used by PetscSolver
		// to compute the row and column indices for the Jacobian
		indices[advClusterIdx] = currId;

		// Compute the partial derivatives for advection of this cluster as
		// explained in the description of this method
		val[advClusterIdx * 2] = -(3.0 * sinkStrength * diffCoeff) /
			(kBoltzmann * cluster.getTemperature(ix + 1) * hxRight *
				pow(pos[0] - location, 4)) *
			advectionGrid[iz + 1][iy + 1][ix + 1][advClusterIdx]; // middle
		val[advClusterIdx * 2] += (3.0 * sinkStrength) *
			(cluster.getDiffusionCoefficient(ix + 2) /
					cluster.getTemperature(ix + 2) -
				cluster.getDiffusionCoefficient(ix + 1) /
					cluster.getTemperature(ix + 1)) /
			(kBoltzmann * hxRight * pow(pos[0] - location, 4)) *
			advectionGrid[iz + 1][iy + 1][ix + 1][advClusterIdx]; // middle
		val[(advClusterIdx * 2) + 1] = (3.0 * sinkStrength * diffCoeff) /
			(kBoltzmann * cluster.getTemperature(ix + 1) * hxRight *
				pow(pos[0] - location + hxRight, 4)) *
			advectionGrid[iz + 1][iy + 1][ix + 2][advClusterIdx]; // right

		++advClusterIdx;
	}

	return;
}

} /* end namespace advection */
} /* end namespace core */
} /* end namespace xolotl */
