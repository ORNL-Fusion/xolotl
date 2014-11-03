// Includes
#include "AdvectionHandler.h"

namespace xolotlCore {

void AdvectionHandler::computeAdvection(std::shared_ptr<PSIClusterReactionNetwork> network,
		double hx, int xi, double *concOffset, double *rightConcOffset,
		double *updatedConcOffset) {
	// Get all the reactant
	auto reactants = network->getAll();
	// Get the number of advecting cluster
	int nAdvec = indexVector.size();

	// Loop on them
	for (int i = 0; i < nAdvec; i++) {
		// Get the advecting cluster
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		// Get the index of the cluster
		int index = cluster->getId() - 1;

		// Get the initial concentrations
		double oldConc = concOffset[index];
		double oldRightConc = rightConcOffset[index];

		// Compute the concentration as explained in the description of the method
		double conc = (3.0 * sinkStrengthVector[i] * cluster->getDiffusionCoefficient())
				* ((oldRightConc / pow((xi + 1) * hx, 4)) - (oldConc / pow(xi * hx, 4)))
				/ (xolotlCore::kBoltzmann * cluster->getTemperature() * hx);

//		double conc = sinkStrengthVector[i] * (oldRightConc - oldConc) / hx;

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;
	}

	return;
}

void AdvectionHandler::computePartialsForAdvection(
		std::shared_ptr<PSIClusterReactionNetwork> network,
		double hx, double *val, int *row, int *col, int xi,
		int xs) {
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
		// Get the sink strenght value
		double sinkStrength = sinkStrengthVector[i];

		// Set the row and column indices. These indices are computed
		// by using xi and xi-1, and the arrays are shifted to
		// (xs+1)*size to properly account for the neighboring ghost
		// cells.

		// Set the row index
		row[i] = (xi - xs + 1) * size + index;

		// Set the columns indices
		col[i * 2] = ((xi - 1) - xs + 1) * size + index;
		col[(i * 2) + 1] = (xi - xs + 1) * size + index;

		// Compute the partial derivatives for advection of this cluster as
		// explained in the description of this method
		val[i * 2] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* hx * pow(xi * hx, 4));
		val[(i * 2) + 1] = -(3.0 * sinkStrength * diffCoeff)
								/ (xolotlCore::kBoltzmann * cluster->getTemperature()
										* hx * pow(xi * hx, 4));

//		val[i * 2] = sinkStrength / hx;
//		val[(i * 2) + 1] = - sinkStrength / hx;
	}

	return;
}

}/* end namespace xolotlCore */
