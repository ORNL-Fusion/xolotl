// Includes
#include "Diffusion2DHandler.h"

namespace xolotlCore {

void Diffusion2DHandler::computeDiffusion(PSIClusterReactionNetwork *network,
		double **concVector, double *updatedConcOffset,
		double sx, double sy, double sz) {
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
		double oldConc = concVector[0][index]; // middle
		double oldLeftConc = concVector[1][index]; // left
		double oldRightConc = concVector[2][index]; // right
		double oldBottomConc = concVector[3][index]; // bottom
		double oldTopConc = concVector[4][index]; // top

		// Use a simple midpoint stencil to compute the concentration
		double conc = cluster->getDiffusionCoefficient()
				* (sx * (oldLeftConc + oldRightConc - 2.0 * oldConc)
						+ sy * (oldBottomConc + oldTopConc - 2.0 * oldConc));

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;
	}

	return;
}

void Diffusion2DHandler::computePartialsForDiffusion(
		PSIClusterReactionNetwork *network,
		double *val, int *indices,
		double sx, double sy, double sz) {
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
		val[i * 5] = -2.0 * diffCoeff * (sx + sy); // middle
		val[(i * 5) + 1] = diffCoeff * sx; // left
		val[(i * 5) + 2] = diffCoeff * sx; // right
		val[(i * 5) + 3] = diffCoeff * sy; // bottom
		val[(i * 5) + 4] = diffCoeff * sy; // top
	}

	return;
}

}/* end namespace xolotlCore */
