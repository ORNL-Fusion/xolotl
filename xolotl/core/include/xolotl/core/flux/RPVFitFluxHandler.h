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

		// Set the flux index corresponding the single void cluster here
		NetworkType::Composition comp = NetworkType::Composition::zero();
		comp[NetworkType::Species::V] = 1;
		auto cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe single void cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::V] = 2;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe V2 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::V] = 3;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe V3 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::V] = 4;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe V4 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::V] = 5;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe V5 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::V] = 9;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe V9 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		// Interstitials
		comp[NetworkType::Species::V] = 0;
		comp[NetworkType::Species::I] = 1;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe single interstitial cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::I] = 2;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe I2 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::I] = 3;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe I3 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::I] = 4;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe I4 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::I] = 5;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe I5 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::I] = 6;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe I6 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::I] = 7;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe I7 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::I] = 8;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe I8 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::I] = 9;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe I9 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::I] = 12;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe I12 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::I] = 16;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe I16 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::I] = 20;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe I20 cluster is not "
									 "present in the network, "
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
			fluxIndices.data(), 18};
		Kokkos::View<IdType*> ids{"Flux Indices", 18};
		deep_copy(ids, ids_h);
		Kokkos::parallel_for(
			1, KOKKOS_LAMBDA(std::size_t) {
				updatedConcOffset[ids[0]] += scalingFactor * 11.776939094; // V1
				updatedConcOffset[ids[1]] += scalingFactor * 1.815018938; // V2
				updatedConcOffset[ids[2]] +=
					scalingFactor * 0.31463055266666667; // V3
				updatedConcOffset[ids[3]] += scalingFactor * 0.191413343; // V4
				updatedConcOffset[ids[4]] += scalingFactor * 0.07611408; // V5
				updatedConcOffset[ids[5]] += scalingFactor * 0.0382584; // V9
				updatedConcOffset[ids[6]] += scalingFactor * 5.940341464; // I1
				updatedConcOffset[ids[7]] += scalingFactor * 2.186892002; // I2
				updatedConcOffset[ids[8]] +=
					scalingFactor * 0.70030164666666661; // I3
				updatedConcOffset[ids[9]] += scalingFactor * 0.447750235; // I4
				updatedConcOffset[ids[10]] += scalingFactor * 0.261860134; // I5
				updatedConcOffset[ids[11]] += scalingFactor * 0.16145861; // I6
				updatedConcOffset[ids[12]] +=
					scalingFactor * 0.094966872285714293; // I7
				updatedConcOffset[ids[13]] += scalingFactor * 0.02310606; // I8
				updatedConcOffset[ids[14]] +=
					scalingFactor * 0.025807277555555556; // I9
				updatedConcOffset[ids[15]] +=
					scalingFactor * 0.0096645931666666674; // I12
				updatedConcOffset[ids[16]] += scalingFactor * 0.00611631; // I16
				updatedConcOffset[ids[17]] +=
					scalingFactor * 0.0030808080; // I20
			});
	}
};
// end class RPVFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl
