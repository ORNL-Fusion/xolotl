#pragma once

#include <cmath>

#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/FeReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the IFluxHandler interface to calculate the incident
 * fluxes for an iron material.
 */
class FeFitFluxHandler : public FluxHandler
{
private:
	/**
	 * \see FluxHandler.h
	 */
	double
	FitFunction(double x)
	{
		// Not actually used
		return 0.0;
	}

public:
	/**
	 * The constructor
	 */
	FeFitFluxHandler(const options::IOptions& options) : FluxHandler(options)
	{
	}

	/**
	 * The Destructor
	 */
	~FeFitFluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid)
	{
		// Only defined in 0D
		if (xGrid.size() == 0) {
			// Add an empty vector
			std::vector<double> tempVector;
			incidentFluxVec.push_back(tempVector);
		}

		using NetworkType = network::FeReactionNetwork;
		auto feNetwork = dynamic_cast<NetworkType*>(&network);

		// Set the flux index corresponding the the single helium cluster here
		NetworkType::Composition comp = NetworkType::Composition::zero();
		comp[NetworkType::Species::He] = 1;
		auto cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe single helium cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		// Look for interstitial now
		comp[NetworkType::Species::He] = 0;
		comp[NetworkType::Species::I] = 1;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe single interstitial cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		// Look for vacancies now
		comp[NetworkType::Species::I] = 0;
		comp[NetworkType::Species::V] = 1;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe single vacancy cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		comp[NetworkType::Species::V] = 2;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe double vacancy cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		comp[NetworkType::Species::V] = 3;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe triple vacancy cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		comp[NetworkType::Species::V] = 4;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe quadruple vacancy cluster is not present "
				"in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		comp[NetworkType::Species::V] = 5;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nVacancy 5 cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		comp[NetworkType::Species::V] = 9;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nVacancy 9 cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		return;
	}

	/**
	 * \see IFluxHandler.h
	 */
    ////////////////////////////////////////////////////////////////////////////
    // DELETEME
	void
	computeIncidentFlux(
		double currentTime, double* updatedConcOffset, int xi, int surfacePos)
	{
		// Define only for a 0D case
		if (incidentFluxVec[0].size() == 0) {
			updatedConcOffset[fluxIndices[0]] += 2.11e-11; // He1
			updatedConcOffset[fluxIndices[1]] += 1.49e-05; // I1
			updatedConcOffset[fluxIndices[2]] += 9.91e-06; // V1
			updatedConcOffset[fluxIndices[3]] += 1.51e-06; // V2
			updatedConcOffset[fluxIndices[4]] += 2.60e-07; // V3
			updatedConcOffset[fluxIndices[5]] += 1.58e-07; // V4
			updatedConcOffset[fluxIndices[6]] += 6.29e-08; // V5
			updatedConcOffset[fluxIndices[7]] += 3.16e-08; // V9
		}

		else {
			throw std::runtime_error(
				"\nThe iron problem is not defined for more than 0D!");
		}

		return;
	}
    ////////////////////////////////////////////////////////////////////////////

	void
	computeIncidentFlux(double currentTime,
		Kokkos::View<double*> updatedConcOffset, int xi,
		int surfacePos) override
	{
		// Define only for a 0D case
		if (incidentFluxVec[0].size() != 0) {
			throw std::runtime_error(
				"\nThe iron problem is not defined for more than 0D!");
		}

		Kokkos::View<IdType*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> ids_h{
			fluxIndices.data(), 8};
		Kokkos::View<IdType*> ids{"Flux Indices", 8};
		deep_copy(ids, ids_h);
		Kokkos::parallel_for(
			1, KOKKOS_LAMBDA(std::size_t) {
				updatedConcOffset[ids[0]] += 2.11e-11; // He1
				updatedConcOffset[ids[1]] += 1.49e-05; // I1
				updatedConcOffset[ids[2]] += 9.91e-06; // V1
				updatedConcOffset[ids[3]] += 1.51e-06; // V2
				updatedConcOffset[ids[4]] += 2.60e-07; // V3
				updatedConcOffset[ids[5]] += 1.58e-07; // V4
				updatedConcOffset[ids[6]] += 6.29e-08; // V5
				updatedConcOffset[ids[7]] += 3.16e-08; // V9
			});
	}
};
// end class FeFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl
