// Includes
#include "XGBAdvectionHandler.h"

namespace xolotlCore {

void XGBAdvectionHandler::initialize(IReactionNetwork *network,
		int *ofill) {
	// Get all the reactants and their number
	auto reactants = network->getAll();
	int networkSize = reactants->size();
	int dof = network->getDOF();

	// Clear the index and sink strength vectors
	indexVector.clear();
	sinkStrengthVector.clear();

	// Loop on all the reactants
	for (int i = 0; i < networkSize; i++) {
		// Get the i-th cluster
		auto cluster = (PSICluster *) reactants->at(i);
		// Get its diffusion coefficient
		double diffFactor = cluster->getDiffusionFactor();

		// Don't do anything if the diffusion factor is 0.0
		if (xolotlCore::equal(diffFactor, 0.0))
			continue;

		// Keep only the helium clusters
		if (cluster->getType() != heType)
			continue;

		// Get its size
		int heSize = cluster->getSize();

		// Switch on the size to get the sink strength (in eV.nm3)
		double sinkStrength = 0.0;
		switch (heSize) {
		case 1:
			sinkStrength = 0.54e-3;
			break;
		case 2:
			sinkStrength = 1.01e-3;
			break;
		case 3:
			sinkStrength = 3.03e-3;
			break;
		case 4:
			sinkStrength = 3.93e-3;
			break;
		case 5:
			sinkStrength = 7.24e-3;
			break;
		case 6:
			sinkStrength = 10.82e-3;
			break;
		case 7:
			sinkStrength = 19.26e-3;
			break;
		}

		// If the sink strength is still 0.0, this cluster is not advecting
		if (xolotlCore::equal(sinkStrength, 0.0))
			continue;

		// Add its index (i) to the vector of indices
		indexVector.push_back(i);

		// Add the sink strength to the vector
		sinkStrengthVector.push_back(sinkStrength);

		// Set the off-diagonal part for the Jacobian to 1
		// Get its id
		int index = cluster->getId() - 1;
		// Set the ofill value to 1 for this cluster
		ofill[index * dof + index] = 1;
	}

	return;
}

void XGBAdvectionHandler::computeAdvection(IReactionNetwork *network,
		std::vector<double> &pos, double **concVector,
		double *updatedConcOffset, double hxLeft, double hxRight, int ix,
		double hy, int iy, double hz, int iz) {
	// Get all the reactant
	auto reactants = network->getAll();
	// Get the number of advecting cluster
	int nAdvec = indexVector.size();

	// Loop on the advecting clusters
	for (int i = 0; i < nAdvec; i++) {
		// Get a specific one and its index
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		int index = cluster->getId() - 1;

		// If we are on the sink, the behavior is not the same
		// Both sides are giving their concentrations to the center
		if (isPointOnSink(pos)) {
			double oldLeftConc = concVector[1][index]; // left
			double oldRightConc = concVector[2][index]; // right

			double conc = (3.0 * sinkStrengthVector[i]
					* cluster->getDiffusionCoefficient())
					* ((oldLeftConc / pow(hxLeft, 5))
							+ (oldRightConc / pow(hxRight, 5)))
					/ (xolotlCore::kBoltzmann * cluster->getTemperature());

			// Update the concentration of the cluster
			updatedConcOffset[index] += conc;
		}
		// Here we are NOT on the GB sink
		else {
			// Get the initial concentrations
			double oldConc = concVector[0][index]; // middle
			double oldRightConc = concVector[2 * (pos[0] > location)
					+ 1 * (pos[0] < location)][index]; // left or right

			// Get the a=d and b=d+h positions
			double a = fabs(location - pos[0]);
			double b = fabs(location - pos[0]) + hxRight * (pos[0] > location)
					+ hxLeft * (pos[0] < location);

			// Compute the concentration as explained in the description of the method
			double conc = (3.0 * sinkStrengthVector[i]
					* cluster->getDiffusionCoefficient())
					* ((oldRightConc / pow(b, 4)) - (oldConc / pow(a, 4)))
					/ (xolotlCore::kBoltzmann * cluster->getTemperature()
							* (hxRight * (pos[0] > location)
									+ hxLeft * (pos[0] < location)));

			// Update the concentration of the cluster
			updatedConcOffset[index] += conc;
		}
	}

	return;
}

void XGBAdvectionHandler::computePartialsForAdvection(
		IReactionNetwork *network, double *val, int *indices,
		std::vector<double> &pos, double hxLeft, double hxRight, int ix,
		double hy, int iy, double hz, int iz) {
	// Get all the reactant
	auto reactants = network->getAll();
	// Get the number of advecting cluster
	int nAdvec = indexVector.size();

	// Loop on the advecting clusters
	for (int i = 0; i < nAdvec; i++) {
		// Get a specific one and its index
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		int index = cluster->getId() - 1;
		// Get the diffusion coefficient of the cluster
		double diffCoeff = cluster->getDiffusionCoefficient();
		// Get the sink strength value
		double sinkStrength = sinkStrengthVector[i];

		// Set the cluster index that will be used by PetscSolver
		// to compute the row and column indices for the Jacobian
		indices[i] = index;

		// If we are on the sink, the partial derivatives are not the same
		// Both sides are giving their concentrations to the center
		if (isPointOnSink(pos)) {
			// 1D case
			if (dimension == 1) {
				val[i * 2] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* pow(hxLeft, 5)); // left
				val[(i * 2) + 1] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* pow(hxRight, 5)); // right
			}
		}
		// Here we are NOT on the GB sink
		else {
			// Get the a=d and b=d+h positions
			double a = fabs(location - pos[0]);
			double b = fabs(location - pos[0]) + hxRight * (pos[0] > location)
					+ hxLeft * (pos[0] < location);

			// Compute the partial derivatives for advection of this cluster as
			// explained in the description of this method
			val[i * 2] = -(3.0 * sinkStrength * diffCoeff)
					/ (xolotlCore::kBoltzmann * cluster->getTemperature()
							* pow(a, 4)
							* (hxRight * (pos[0] > location)
									+ hxLeft * (pos[0] < location))); // middle
			val[(i * 2) + 1] = (3.0 * sinkStrength * diffCoeff)
					/ (xolotlCore::kBoltzmann * cluster->getTemperature()
							* pow(b, 4)
							* (hxRight * (pos[0] > location)
									+ hxLeft * (pos[0] < location))); // left or right
		}
	}

	return;
}

std::vector<int> XGBAdvectionHandler::getStencilForAdvection(
		std::vector<double> &pos) {
	// The first index is positive by convention if we are on the sink
	if (isPointOnSink(pos))
		return {1, 0, 0};
	// The first index is positive if pos[0] > location
	// negative if pos[0] < location
	return {(pos[0] > location) - (pos[0] < location), 0, 0};
}

bool XGBAdvectionHandler::isPointOnSink(std::vector<double> &pos) {
	// Return true if pos[0] is equal to location
	return fabs(location - pos[0]) < 0.001;
}

}/* end namespace xolotlCore */
