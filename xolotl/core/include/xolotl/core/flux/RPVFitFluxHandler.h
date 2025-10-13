#pragma once

#include <cmath>

#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/RPVReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the IFluxHandler interface to calculate the incident
 * fluxes for an Fe BCC material.
 */
class RPVFitFluxHandler : public FluxHandler
{
private:
	/**
	 * \see FluxHandler.h
	 */
	double
	FitFunction(double x) override
	{
		// Not actually used
		return 0.0;
	}

	/**
	 * Scaling the generation term
	 */
	double scalingFactor = 1.0;

public:
	/**
	 * The constructor
	 */
	RPVFitFluxHandler(const options::IOptions& options) : FluxHandler(options)
	{
	}

	/**
	 * The Destructor
	 */
	~RPVFitFluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid) override
	{
		// Only defined in 0D
		if (xGrid.size() == 0) {
			// Add an empty vector
			std::vector<double> tempVector;
			incidentFluxVec.push_back(tempVector);
		}

		// To scale the values to the default dose rate of 1.0e-6 dpa s-1
		scalingFactor = getFluxAmplitude() * 1.0e6;

		using NetworkType = network::RPVReactionNetwork;
		auto feNetwork = dynamic_cast<NetworkType*>(&network);

		// Set the flux index corresponding the the single helium cluster here
		NetworkType::Composition comp = NetworkType::Composition::zero();
		comp[NetworkType::Species::I] = 1;
		auto cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
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
	void
	computeIncidentFlux(double currentTime, Kokkos::View<const double*>,
		Kokkos::View<double*> updatedConcOffset, int xi,
		int surfacePos) override
	{
		// Define only for a 0D case
		if (incidentFluxVec[0].size() != 0) {
			throw std::runtime_error(
				"\nThe iron problem is not defined for more than 0D!");
		}

		Kokkos::View<IdType*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> ids_h{
			fluxIndices.data(), 7};
		Kokkos::View<IdType*> ids{"Flux Indices", 7};
		deep_copy(ids, ids_h);
		Kokkos::parallel_for(
			1, KOKKOS_LAMBDA(std::size_t) {
				updatedConcOffset[ids[0]] += 1.49e-05 * scalingFactor; // I1
				updatedConcOffset[ids[1]] += 9.91e-06 * scalingFactor; // V1
				updatedConcOffset[ids[2]] += 1.51e-06 * scalingFactor; // V2
				updatedConcOffset[ids[3]] += 2.60e-07 * scalingFactor; // V3
				updatedConcOffset[ids[4]] += 1.58e-07 * scalingFactor; // V4
				updatedConcOffset[ids[5]] += 6.29e-08 * scalingFactor; // V5
				updatedConcOffset[ids[6]] += 3.16e-08 * scalingFactor; // V9
			});
	}
};
// end class RPVFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl
