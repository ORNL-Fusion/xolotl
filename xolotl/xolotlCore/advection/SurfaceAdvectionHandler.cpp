// Includes
#include "SurfaceAdvectionHandler.h"

namespace xolotlCore {

void SurfaceAdvectionHandler::initializeAdvectionGrid(
		std::vector<IAdvectionHandler *> advectionHandlers,
		std::vector<double> grid, int ny, double hy, int nz, double hz) {

	// Get the number of advecting clusters
	int nAdvec = advectingClusters.size();

	// Get the size of the grid in the depth direction
	int nx = grid.size();

	// Initialize the diffusion grid with true everywhere
	advectionGrid.clear();
	// Initialize it to True
	for (int k = 0; k < nz + 2; k++) {
		std::vector<std::vector<std::vector<bool> > > tempGridTer;
		for (int j = 0; j < ny + 2; j++) {
			std::vector<std::vector<bool> > tempGridBis;
			for (int i = 0; i < nx + 2; i++) {
				tempGridBis.emplace_back(nAdvec, true);
			}
			tempGridTer.push_back(tempGridBis);
		}
		advectionGrid.push_back(tempGridTer);
	}

	// Initialize the grid position
	Point3D gridPosition { 0.0, 0.0, 0.0 };

	// Consider each advection handler.
	for (auto const& currAdvecHandler : advectionHandlers) {

		// Get the list of advecting clusters
		auto const& otherAdvecClusters =
				currAdvecHandler->getAdvectingClusters();

		// Loop on the spatial grid
		for (int k = -1; k < nz + 1; k++) {
			// Set the grid position
			gridPosition[2] = hz * (double) k;
			for (int j = -1; j < ny + 1; j++) {
				// Set the grid position
				gridPosition[1] = hy * (double) j;
				for (int i = -1; i < nx + 1; i++) {
					// Set the grid position
					if (i == -1)
						gridPosition[0] = -1.0;
					else if (i == nx)
						gridPosition[0] = grid[i - 1] + 1.0;
					else
						gridPosition[0] = grid[i];

					// Check if we are on a sink
					if (currAdvecHandler->isPointOnSink(gridPosition)) {
						// We have to find the corresponding index in the diffusion
						// index vector
						for (int m = 0; m < otherAdvecClusters.size(); m++) {
							// Initialize n the index in the diffusion index vector
							// TODO can we do this with std::find or std::find_if?
							int n = 0;
							while (n < nAdvec) {
								if (&(advectingClusters[n])
										== &(otherAdvecClusters[m])) {
									break;
								}
								n++;
							}
							// Set this diffusion grid value to false
							advectionGrid[k + 1][j + 1][i + 1][n] = false;
						}
					}
				}
			}
		}
	}

	return;
}

void SurfaceAdvectionHandler::computeAdvection(const IReactionNetwork& network,
		const Point3D& pos, double **concVector, double *updatedConcOffset,
		double hxLeft, double hxRight, int ix, double hy, int iy, double hz,
		int iz) const {

	// Get the number of advecting cluster
	int nAdvec = advectingClusters.size();

	// Consider each advecting cluster
	// TODO Maintaining a separate index assumes that advectingClusters is
	// visited in same order as advectionGrid array for given point
	// and the sinkStrengthVector.
	// Currently true with C++11, but we'd like to be able to visit the
	// advecting clusters in any order (so that we can parallelize).
	// Maybe with a zip? or a std::transform?
	int advClusterIdx = 0;
	for (IReactant const& currReactant : advectingClusters) {
		// Get a specific one and its index
		auto const& cluster = static_cast<PSICluster const&>(currReactant);
		int index = cluster.getId() - 1;

		// Get the initial concentrations
		double oldConc = concVector[0][index]
				* advectionGrid[iz + 1][iy + 1][ix + 1][advClusterIdx]; // middle
		double oldRightConc = concVector[2][index]
				* advectionGrid[iz + 1][iy + 1][ix + 2][advClusterIdx]; // right

		// Compute the concentration as explained in the description of the method
		double conc = (3.0 * sinkStrengthVector[advClusterIdx]
				* cluster.getDiffusionCoefficient())
				* ((oldRightConc / pow(pos[0] - location + hxRight, 4))
						- (oldConc / pow(pos[0] - location, 4)))
				/ (xolotlCore::kBoltzmann * cluster.getTemperature() * hxRight);

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;

		++advClusterIdx;
	}

	return;
}

void SurfaceAdvectionHandler::computePartialsForAdvection(
		const IReactionNetwork& network, double *val, int *indices,
		const Point3D& pos, double hxLeft, double hxRight, int ix, double hy,
		int iy, double hz, int iz) const {

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
	for (IReactant const& currReactant : advectingClusters) {
		// Get a specific one and its index
		auto const& cluster = static_cast<PSICluster const&>(currReactant);
		int index = cluster.getId() - 1;
		// Get the diffusion coefficient of the cluster
		double diffCoeff = cluster.getDiffusionCoefficient();
		// Get the sink strength value
		double sinkStrength = sinkStrengthVector[advClusterIdx];

		// Set the cluster index that will be used by PetscSolver
		// to compute the row and column indices for the Jacobian
		indices[advClusterIdx] = index;

		// Compute the partial derivatives for advection of this cluster as
		// explained in the description of this method
		val[advClusterIdx * 2] = -(3.0 * sinkStrength * diffCoeff)
				/ (xolotlCore::kBoltzmann * cluster.getTemperature() * hxRight
						* pow(pos[0] - location, 4))
				* advectionGrid[iz + 1][iy + 1][ix + 1][advClusterIdx]; // middle
		val[(advClusterIdx * 2) + 1] = (3.0 * sinkStrength * diffCoeff)
				/ (xolotlCore::kBoltzmann * cluster.getTemperature() * hxRight
						* pow(pos[0] - location + hxRight, 4))
				* advectionGrid[iz + 1][iy + 1][ix + 2][advClusterIdx]; // right

		++advClusterIdx;
	}

	return;
}

}/* end namespace xolotlCore */
