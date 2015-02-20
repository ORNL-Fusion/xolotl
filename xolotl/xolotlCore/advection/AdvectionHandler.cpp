// Includes
#include "AdvectionHandler.h"

namespace xolotlCore {

void AdvectionHandler::computeAdvection(PSIClusterReactionNetwork *network,
		double hx, std::vector<double> &pos, double **concVector,
		double *updatedConcOffset) {
	// Get all the reactant
	auto reactants = network->getAll();
	// Get the number of advecting cluster
	int nAdvec = indexVector.size();

	// Get the number of degrees of freedom which is the size of the network
	int dof = reactants->size();

	// Loop on them
	for (int i = 0; i < nAdvec; i++) {
		// Get the advecting cluster
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		// Get the index of the cluster
		int index = cluster->getId() - 1;

		// Get the initial concentrations
		double oldConc = concVector[0][index]; // middle
		double oldRightConc = concVector[2][index]; // right

		// Compute the concentration as explained in the description of the method
		double conc = (3.0 * sinkStrengthVector[i] * cluster->getDiffusionCoefficient())
				* ((oldRightConc / pow(pos[0] + hx, 4)) - (oldConc / pow(pos[0], 4)))
				/ (xolotlCore::kBoltzmann * cluster->getTemperature() * hx);

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;
	}

	return;
}

void AdvectionHandler::computePartialsForAdvection(
		PSIClusterReactionNetwork *network,
		double hx, double *val, int *indices, std::vector<double> &pos) {
	// Get all the reactant
	auto reactants = network->getAll();
	// And the size of the network
	int size = reactants->size();
	// Get the number of diffusing cluster
	int nAdvec = indexVector.size();

	// Loop on them
	for (int i = 0; i < nAdvec; i++) {
		// Get the diffusing cluster
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		// Get the index of the cluster
		int index = cluster->getId() - 1;
		// Get the diffusion coefficient of the cluster
		double diffCoeff = cluster->getDiffusionCoefficient();
		// Get the sink strength value
		double sinkStrength = sinkStrengthVector[i];

		// Set the cluster index that will be used by PetscSolver
		// to compute the row and column indices for the Jacobian
		indices[i] = index;

		// Compute the partial derivatives for advection of this cluster as
		// explained in the description of this method
		val[i * 2] = -(3.0 * sinkStrength * diffCoeff)
				/ (xolotlCore::kBoltzmann * cluster->getTemperature()
						* hx * pow(pos[0], 4)); // middle
		val[(i * 2) + 1] = (3.0 * sinkStrength * diffCoeff)
				/ (xolotlCore::kBoltzmann * cluster->getTemperature()
						* hx * pow(pos[0] + hx, 4)); // right
	}

	return;
}

}/* end namespace xolotlCore */
