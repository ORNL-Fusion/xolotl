// Includes
#include "ZGBAdvectionHandler.h"

namespace xolotlCore {

void ZGBAdvectionHandler::initialize(const IReactionNetwork& network,
                    IReactionNetwork::SparseFillMap& ofillMap) {

	int dof = network.getDOF();

	// Clear the index and sink strength vectors
	advectingClusters.clear();
	sinkStrengthVector.clear();

	// Loop on all the reactants
	for (IReactant const& currReactant : network.getAll()) {

		auto const& cluster = static_cast<PSICluster const&>(currReactant);
		// Get its diffusion coefficient
		double diffFactor = cluster.getDiffusionFactor();

		// Don't do anything if the diffusion factor is 0.0
		if (xolotlCore::equal(diffFactor, 0.0))
			continue;

		// Keep only the helium clusters
		if (cluster.getType() != ReactantType::He)
			continue;

		// Get its size
		int heSize = cluster.getSize();

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

		// Note that the current cluster is advecting.
		advectingClusters.emplace_back(cluster);

		// Add the sink strength to the vector
		sinkStrengthVector.push_back(sinkStrength);

		// Set the off-diagonal part for the Jacobian to 1
		// Get its id
		int index = cluster.getId() - 1;
		// Set the ofill value to 1 for this cluster
        ofillMap[index].emplace_back(index);
	}

	return;
}

void ZGBAdvectionHandler::computeAdvection(const IReactionNetwork& network,
		const Point3D& pos, double **concVector, double *updatedConcOffset,
		double hxLeft, double hxRight, int ix, double hy, int iy, double hz,
		int iz) const {

	// Consider each advecting cluster.
	// TODO Maintaining a separate index assumes that advectingClusters is
	// visited in same order as advectionGrid array for given point
	// and the sinkStrengthVector.
	// Currently true with C++11, but we'd like to be able to visit the
	// advecting clusters in any order (so that we can parallelize).
	// Maybe with a zip? or a std::transform?
	int advClusterIdx = 0;
	for (IReactant const& currReactant : advectingClusters) {

		auto const& cluster = static_cast<PSICluster const&>(currReactant);
		int index = cluster.getId() - 1;

		// If we are on the sink, the behavior is not the same
		// Both sides are giving their concentrations to the center
		if (isPointOnSink(pos)) {
			double oldFrontConc = concVector[5][index]; // front
			double oldBackConc = concVector[6][index]; // back

			double conc = (3.0 * sinkStrengthVector[advClusterIdx]
					* cluster.getDiffusionCoefficient())
					* ((oldFrontConc + oldBackConc) / pow(hz, 5))
					/ (xolotlCore::kBoltzmann * cluster.getTemperature());

			// Update the concentration of the cluster
			updatedConcOffset[index] += conc;
		}
		// Here we are NOT on the GB sink
		else {
			// Get the initial concentrations
			double oldConc = concVector[0][index]; // middle
			double oldRightConc = concVector[6 * (pos[2] > location)
					+ 5 * (pos[2] < location)][index]; // back or front

			// Get the a=d and b=d+h positions
			double a = fabs(location - pos[2]);
			double b = fabs(location - pos[2]) + hz;

			// Compute the concentration as explained in the description of the method
			double conc = (3.0 * sinkStrengthVector[advClusterIdx]
					* cluster.getDiffusionCoefficient())
					* ((oldRightConc / pow(b, 4)) - (oldConc / pow(a, 4)))
					/ (xolotlCore::kBoltzmann * cluster.getTemperature() * hz);

			// Update the concentration of the cluster
			updatedConcOffset[index] += conc;
		}

		++advClusterIdx;
	}

	return;
}

void ZGBAdvectionHandler::computePartialsForAdvection(
		const IReactionNetwork& network, double *val, int *indices,
		const Point3D& pos, double hxLeft, double hxRight, int ix, double hy,
		int iy, double hz, int iz) const {

	// Loop on the advecting clusters
	// TODO Maintaining a separate index assumes that advectingClusters is
	// visited in same order as advectionGrid array for given point
	// and the sinkStrengthVector.
	// Currently true with C++11, but we'd like to be able to visit the
	// advecting clusters in any order (so that we can parallelize).
	// Maybe with a zip? or a std::transform?
	int advClusterIdx = 0;
	for (IReactant const& currReactant : advectingClusters) {

		auto const& cluster = static_cast<PSICluster const&>(currReactant);

		int index = cluster.getId() - 1;
		// Get the diffusion coefficient of the cluster
		double diffCoeff = cluster.getDiffusionCoefficient();
		// Get the sink strength value
		double sinkStrength = sinkStrengthVector[advClusterIdx];

		// Set the cluster index that will be used by PetscSolver
		// to compute the row and column indices for the Jacobian
		indices[advClusterIdx] = index;

		// If we are on the sink, the partial derivatives are not the same
		// Both sides are giving their concentrations to the center
		if (isPointOnSink(pos)) {
			val[advClusterIdx * 2] = (3.0 * sinkStrength * diffCoeff)
					/ (xolotlCore::kBoltzmann * cluster.getTemperature()
							* pow(hz, 5)); // back or front
			val[(advClusterIdx * 2) + 1] = val[advClusterIdx * 2]; // back or front
		}
		// Here we are NOT on the GB sink
		else {
			// Get the a=d and b=d+h positions
			double a = fabs(location - pos[2]);
			double b = fabs(location - pos[2]) + hz;

			// Compute the partial derivatives for advection of this cluster as
			// explained in the description of this method
			val[advClusterIdx * 2] = -(3.0 * sinkStrength * diffCoeff)
					/ (xolotlCore::kBoltzmann * cluster.getTemperature() * hz
							* pow(a, 4)); // middle
			val[(advClusterIdx * 2) + 1] = (3.0 * sinkStrength * diffCoeff)
					/ (xolotlCore::kBoltzmann * cluster.getTemperature() * hz
							* pow(b, 4)); // back or front
		}

		++advClusterIdx;
	}

	return;
}

std::array<int, 3> ZGBAdvectionHandler::getStencilForAdvection(
		const Point3D& pos) const {

	// The third index is positive by convention if we are on the sink
	if (isPointOnSink(pos))
		return {0, 0, 1};
	// The third index is positive if pos[2] > location
	// negative if pos[2] < location
	return {0, 0, (pos[2] > location) - (pos[2] < location)};
}

}/* end namespace xolotlCore */
