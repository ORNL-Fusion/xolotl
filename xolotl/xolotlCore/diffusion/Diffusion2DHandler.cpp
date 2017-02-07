// Includes
#include "Diffusion2DHandler.h"

namespace xolotlCore {

void Diffusion2DHandler::initializeDiffusionGrid(std::vector<IAdvectionHandler *> advectionHandlers,
		std::vector<double> grid,
		int ny, double hy, int nz, double hz) {
	// Get the number of diffusing clusters
	int nDiff = indexVector.size();

	// Get the size of the grid in the depth direction
	int nx = grid.size();

	// Initialize the diffusion grid with true everywhere
	diffusionGrid.clear();
	// Initialize it to True
	for (int j = 0; j < ny + 2; j++) {
		std::vector < std::vector<bool> > tempGridBis;
		for (int i = 0; i < nx + 2; i++) {
			std::vector<bool> tempGrid;
			for (int n = 0; n < nDiff; n++) {
				tempGrid.push_back(true);
			}
			tempGridBis.push_back(tempGrid);
		}
		diffusionGrid.push_back(tempGridBis);
	}

	// Initialize the grid position
	std::vector<double> gridPosition = { 0.0, 0.0, 0.0 };

	// Loop on the advection handlers
	for (int l = 0; l < advectionHandlers.size(); l++) {
		// Get the list of advecting clusters
		auto advecVector = advectionHandlers[l]->getIndexVector();

		// Loop on the spatial grid
		for (int j = -1; j < ny + 1; j++) {
			// Set the grid position
			gridPosition[1] = hy * (double) j;
			for (int i = -1; i < nx + 1; i++) {
				// Set the grid position
				if (i == -1) gridPosition[0] = - 1.0;
				else if (i == nx) gridPosition[0] = grid[i-1] + 1.0;
				else gridPosition[0] = grid[i];

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
						diffusionGrid[j+1][i+1][n] = false;
					}
				}
			}
		}
	}

	return;
}

void Diffusion2DHandler::computeDiffusion(IReactionNetwork *network,
		double **concVector, double *updatedConcOffset,
		double hxLeft, double hxRight, int ix,
		double sy, int iy, double, int) {
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
		double oldConc = concVector[0][index] * diffusionGrid[iy+1][ix+1][i]; // middle
		double oldLeftConc = concVector[1][index] * diffusionGrid[iy+1][ix][i]; // left
		double oldRightConc = concVector[2][index] * diffusionGrid[iy+1][ix+2][i]; // right
		double oldBottomConc = concVector[3][index] * diffusionGrid[iy][ix+1][i]; // bottom
		double oldTopConc = concVector[4][index] * diffusionGrid[iy+2][ix+1][i]; // top

		// Use a simple midpoint stencil to compute the concentration
		double conc = cluster->getDiffusionCoefficient()
				* (2.0 * (oldLeftConc + (hxLeft / hxRight) * oldRightConc
								- (1.0 + (hxLeft / hxRight)) * oldConc)
						/ (hxLeft * (hxLeft + hxRight))
						+ sy * (oldBottomConc + oldTopConc - 2.0 * oldConc));

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;
	}

	return;
}

void Diffusion2DHandler::computePartialsForDiffusion(
		IReactionNetwork *network,
		double *val, int *indices, double hxLeft, double hxRight, int ix,
		double sy, int iy, double, int) {
	// Get all the reactant
	auto reactants = network->getAll();
	// Get the number of diffusing cluster
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
		// for the middle, left, right, bottom, and top grid point
		val[i * 5] = - 2.0 * diffCoeff * ((1.0 / (hxLeft * hxRight)) + sy) * diffusionGrid[iy+1][ix+1][i]; // middle
		val[(i * 5) + 1] = diffCoeff * 2.0 / (hxLeft * (hxLeft + hxRight)) * diffusionGrid[iy+1][ix][i]; // left
		val[(i * 5) + 2] = diffCoeff * 2.0 / (hxRight * (hxLeft + hxRight)) * diffusionGrid[iy+1][ix+2][i]; // right
		val[(i * 5) + 3] = diffCoeff * sy * diffusionGrid[iy][ix+1][i]; // bottom
		val[(i * 5) + 4] = diffCoeff * sy * diffusionGrid[iy+2][ix+1][i]; // top
	}

	return;
}

}/* end namespace xolotlCore */
