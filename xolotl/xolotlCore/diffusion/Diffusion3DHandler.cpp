// Includes
#include "Diffusion3DHandler.h"

namespace xolotlCore {

void Diffusion3DHandler::initializeDiffusionGrid(std::vector<IAdvectionHandler *> advectionHandlers,
		std::vector<double> grid,
		int ny, double hy, int nz, double hz) {
	// Get the number of diffusing clusters
	int nDiff = indexVector.size();

	// Get the size of the grid in the depth direction
	int nx = grid.size();

	// Initialize the diffusion grid with true everywhere
	diffusionGrid.clear();
	// Initialize it to True
	for (int k = 0; k < nz + 2; k++) {
		std::vector < std::vector < std::vector<bool> > > tempGridTer;
		for (int j = 0; j < ny + 2; j++) {
			std::vector < std::vector<bool> > tempGridBis;
			for (int i = 0; i < nx; i++) {
				std::vector<bool> tempGrid;
				for (int n = 0; n < nDiff; n++) {
					tempGrid.push_back(true);
				}
				tempGridBis.push_back(tempGrid);
			}
			tempGridTer.push_back(tempGridBis);
		}
		diffusionGrid.push_back(tempGridTer);
	}

	// Initialize the grid position
	std::vector<double> gridPosition = { 0.0, 0.0, 0.0 };

	// Loop on the advection handlers
	for (int l = 0; l < advectionHandlers.size(); l++) {
		// Get the list of advecting clusters
		auto advecVector = advectionHandlers[l]->getIndexVector();

		// Loop on the spatial grid
		for (int k = -1; k < nz + 1; k++) {
			// Set the grid position
			gridPosition[2] = hz * (double) k;
			for (int j = -1; j < ny + 1; j++) {
				// Set the grid position
				gridPosition[1] = hy * (double) j;
				for (int i = 0; i < nx; i++) {
					// Set the grid position
					gridPosition[0] = grid[i] - grid[1];

					// Check if we are on a sink
					if (advectionHandlers[l]->isPointOnSink(gridPosition)) {
						// We have to find the corresponding index in the diffusion
						// index vector
						for (int m = 0; m < advecVector.size(); m++) {
							// Initialize n the index in the diffusion index vector
							int n = 0;
							while (n < nDiff) {
								if (indexVector[n] == advecVector[m]) break;
								n++;
							}
							// Set this diffusion grid value to false
							diffusionGrid[k+1][j+1][i][n] = false;
						}
					}
				}
			}
		}
	}

	return;
}

void Diffusion3DHandler::computeDiffusion(IReactionNetwork *network,
		double **concVector, double *updatedConcOffset,
		double hxLeft, double hxRight, int ix,
		double sy, int iy, double sz, int iz) {
	// Get all the reactants
	auto reactants = network->getAll();
	// Get the number of diffusing clusters
	int nDiff = indexVector.size();

	// Loop on them
	for (int i = 0; i < nDiff; i++) {
		// Get the diffusing cluster and its index
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		int index = cluster->getId() - 1;

		// Get the initial concentrations
		double oldConc = concVector[0][index] * diffusionGrid[iz+1][iy+1][ix+1][i]; // middle
		double oldLeftConc = concVector[1][index] * diffusionGrid[iz+1][iy+1][ix][i]; // left
		double oldRightConc = concVector[2][index] * diffusionGrid[iz+1][iy+1][ix+2][i]; // right
		double oldBottomConc = concVector[3][index] * diffusionGrid[iz+1][iy][ix+1][i]; // bottom
		double oldTopConc = concVector[4][index] * diffusionGrid[iz+1][iy+2][ix+1][i]; // top
		double oldFrontConc = concVector[5][index] * diffusionGrid[iz][iy+1][ix+1][i]; // front
		double oldBackConc = concVector[6][index] * diffusionGrid[iz+2][iy+1][ix+1][i]; // back

		// Use a simple midpoint stencil to compute the concentration
		double conc = cluster->getDiffusionCoefficient()
						* (2.0 * (oldLeftConc + (hxLeft / hxRight) * oldRightConc
										- (1.0 + (hxLeft / hxRight)) * oldConc)
								/ (hxLeft * (hxLeft + hxRight))
								+ sy * (oldBottomConc + oldTopConc - 2.0 * oldConc)
								+ sz * (oldFrontConc + oldBackConc - 2.0 * oldConc));

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;
	}

	return;
}

void Diffusion3DHandler::computePartialsForDiffusion(
		IReactionNetwork *network,
		double *val, int *indices, double hxLeft, double hxRight, int ix,
		double sy, int iy, double sz, int iz) {
	// Get all the reactants
	auto reactants = network->getAll();
	// Get the number of diffusing clusters
	int nDiff = indexVector.size();

	// Loop on them
	for (int i = 0; i < nDiff; i++) {
		// Get the diffusing cluster and its index
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		int index = cluster->getId() - 1;
		// Get the diffusion coefficient of the cluster
		double diffCoeff = cluster->getDiffusionCoefficient();

		// Set the cluster index, the PetscSolver will use it to compute
		// the row and column indices for the Jacobian
		indices[i] = index;

		// Compute the partial derivatives for diffusion of this cluster
		// for the middle, left, right, bottom, top, front, and back grid point
		val[i * 7] = - 2.0 * diffCoeff * ((1.0 / (hxLeft * hxRight)) + sy + sz)
				* diffusionGrid[iz+1][iy+1][ix+1][i]; // middle
		val[(i * 7) + 1] = diffCoeff * 2.0 / (hxLeft * (hxLeft + hxRight))
				* diffusionGrid[iz+1][iy+1][ix][i]; // left
		val[(i * 7) + 2] = diffCoeff * 2.0 / (hxRight * (hxLeft + hxRight))
				* diffusionGrid[iz+1][iy+1][ix+2][i]; // right
		val[(i * 7) + 3] = diffCoeff * sy * diffusionGrid[iz+1][iy][ix+1][i]; // bottom
		val[(i * 7) + 4] = diffCoeff * sy * diffusionGrid[iz+1][iy+2][ix+1][i]; // top
		val[(i * 7) + 5] = diffCoeff * sz * diffusionGrid[iz][iy+1][ix+1][i]; // front
		val[(i * 7) + 6] = diffCoeff * sz * diffusionGrid[iz+2][iy+1][ix+1][i]; // back
	}

	return;
}

}/* end namespace xolotlCore */
