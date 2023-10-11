// Includes
#include <xolotl/core/advection/ZGBAdvectionHandler.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace advection
{
void
ZGBAdvectionHandler::initialize(
	network::IReactionNetwork& network, std::vector<RowColPair>& idPairs)
{
	// Clear the index and sink strength vectors
	advectingClusters.clear();
	sinkStrengthVector.clear();

	using NetworkType = network::IPSIReactionNetwork;
	using AmountType = NetworkType::AmountType;

	auto psiNetwork = dynamic_cast<NetworkType*>(&network);
	auto numSpecies = psiNetwork->getSpeciesListSize();
	auto specIdHe = psiNetwork->getHeliumSpeciesId();

	// Initialize the composition
	auto comp = std::vector<AmountType>(numSpecies, 0);

	// Loop on helium clusters from size 1 to 7
	for (std::size_t i = 1; i <= 7; i++) {
		comp[specIdHe()] = i;
		auto clusterId = psiNetwork->findClusterId(comp);

		// Check that the helium cluster is present in the network
		if (clusterId == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe helium cluster of size " +
				std::to_string(i) +
				"is not present in the network, "
				"cannot use the advection option!");
		}

		auto cluster = psiNetwork->getClusterCommon(clusterId);

		// Get its diffusion factor
		double diffFactor = cluster.getDiffusionFactor();

		// Don't do anything if the diffusion factor is 0.0
		if (util::equal(diffFactor, 0.0))
			continue;

		// Switch on the size to get the sink strength (in eV.nm3)
		double sinkStrength = 0.0;
		switch (i) {
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
		if (util::equal(sinkStrength, 0.0))
			continue;

		// Add it to our collection of advecting clusters.
		advectingClusters.emplace_back(clusterId);

		// Add the sink strength to the vector
		sinkStrengthVector.push_back(sinkStrength);

		// Add Jacobian entry for this cluster
		idPairs.push_back({clusterId, clusterId});
	}

	this->syncAdvectingClusters(network);
	this->syncSinkStrengths();
}

void
ZGBAdvectionHandler::computeAdvection(network::IReactionNetwork& network,
	const plsm::SpaceVector<double, 3>& pos, const StencilConcArray& concVector,
	Kokkos::View<double*> updatedConcOffset, double hxLeft, double hxRight,
	int ix, double hy, int iy, double hz, int iz) const
{
	// Consider each advecting cluster.
	// TODO Maintaining a separate index assumes that advectingClusters is
	// visited in same order as advectionGrid array for given point
	// and the sinkStrengthVector.
	// Currently true with C++11, but we'd like to be able to visit the
	// advecting clusters in any order (so that we can parallelize).
	// Maybe with a zip? or a std::transform?

	if (concVector.size() < 7) {
		throw std::runtime_error("Wrong sizse for 3D concentration stencil; "
								 "should be at least 5, got " +
			std::to_string(concVector.size()));
	}
	Kokkos::Array<Kokkos::View<const double*>, 7> concVec = {concVector[0],
		concVector[1], concVector[2], concVector[3], concVector[4],
		concVector[5], concVector[6]};

	auto location_ = location;
	auto clusterIds = this->advClusterIds;
	auto clusters = this->advClusters;
	auto sinkStrengths = this->advSinkStrengths;

	// If we are on the sink, the behavior is not the same
	// Both sides are giving their concentrations to the center
	if (isPointOnSink(pos)) {
		Kokkos::parallel_for(
			clusterIds.size(), KOKKOS_LAMBDA(IdType i) {
				auto currId = clusterIds[i];
				auto cluster = clusters[i];

				double oldFrontConc = concVec[5][currId]; // front
				double oldBackConc = concVec[6][currId]; // back

				double conc = (3.0 * sinkStrengths[i] *
								  cluster.getDiffusionCoefficient(ix + 1)) *
					((oldFrontConc + oldBackConc) / pow(hz, 5)) /
					(kBoltzmann * cluster.getTemperature(ix + 1));

				// Update the concentration of the cluster
				updatedConcOffset[currId] += conc;
			});
	}
	// Here we are NOT on the GB sink
	else {
		Kokkos::parallel_for(
			clusterIds.size(), KOKKOS_LAMBDA(IdType i) {
				auto currId = clusterIds[i];
				auto cluster = clusters[i];

				// Get the initial concentrations
				double oldConc = concVec[0][currId]; // middle
				double oldRightConc = concVec[6 * (pos[2] > location_) +
					5 * (pos[2] < location_)][currId]; // back or front

				// Get the a=d and b=d+h positions
				double a = fabs(location_ - pos[2]);
				double b = fabs(location_ - pos[2]) + hz;

				// Compute the concentration as explained in the description of
				// the method
				double conc = (3.0 * sinkStrengths[i] *
								  cluster.getDiffusionCoefficient(ix + 1)) *
					((oldRightConc / pow(b, 4)) - (oldConc / pow(a, 4))) /
					(kBoltzmann * cluster.getTemperature(ix + 1) * hz);

				// Update the concentration of the cluster
				updatedConcOffset[currId] += conc;
			});
	}
}

void
ZGBAdvectionHandler::computePartialsForAdvection(
	network::IReactionNetwork& network, double* val, IdType* indices,
	const plsm::SpaceVector<double, 3>& pos, double hxLeft, double hxRight,
	int ix, double hy, int iy, double hz, int iz) const
{
	// Loop on the advecting clusters
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

		// If we are on the sink, the partial derivatives are not the same
		// Both sides are giving their concentrations to the center
		if (isPointOnSink(pos)) {
			val[advClusterIdx * 2] = (3.0 * sinkStrength * diffCoeff) /
				(kBoltzmann * cluster.getTemperature(ix + 1) *
					pow(hz, 5)); // back or front
			val[(advClusterIdx * 2) + 1] =
				val[advClusterIdx * 2]; // back or front
		}
		// Here we are NOT on the GB sink
		else {
			// Get the a=d and b=d+h positions
			double a = fabs(location - pos[2]);
			double b = fabs(location - pos[2]) + hz;

			// Compute the partial derivatives for advection of this cluster as
			// explained in the description of this method
			val[advClusterIdx * 2] = -(3.0 * sinkStrength * diffCoeff) /
				(kBoltzmann * cluster.getTemperature(ix + 1) * hz *
					pow(a, 4)); // middle
			val[(advClusterIdx * 2) + 1] = (3.0 * sinkStrength * diffCoeff) /
				(kBoltzmann * cluster.getTemperature(ix + 1) * hz *
					pow(b, 4)); // back or front
		}

		++advClusterIdx;
	}

	return;
}

std::array<int, 3>
ZGBAdvectionHandler::getStencilForAdvection(
	const plsm::SpaceVector<double, 3>& pos) const
{
	// The third index is positive by convention if we are on the sink
	if (isPointOnSink(pos))
		return {0, 0, 1};
	// The third index is positive if pos[2] > location
	// negative if pos[2] < location
	return {0, 0, (pos[2] > location) - (pos[2] < location)};
}

} /* end namespace advection */
} /* end namespace core */
} /* end namespace xolotl */
