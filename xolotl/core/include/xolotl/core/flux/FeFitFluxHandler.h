#pragma once

#include <cmath>

#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/FeReactionNetwork.h>
#include <xolotl/core/network/T91ReactionNetwork.h>

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
	FitFunction(double x) override
	{
		// Not actually used
		return 0.0;
	}

	double hePPM; // ppm / dpa

public:
	/**
	 * The constructor
	 */
	FeFitFluxHandler(const options::IOptions& options) : FluxHandler(options)
	{
		hePPM = options.getHePPM(); // ppm / dpa
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
		std::vector<double> grid) override
	{
		if (auto feNetwork =
				dynamic_cast<network::FeReactionNetwork*>(&network)) {
			addClusters(*feNetwork);
		}
		else if (auto t91Network =
					 dynamic_cast<network::T91ReactionNetwork*>(&network)) {
			addClusters(*t91Network);
		}

		// Copy data to device views
		syncFluxIndices();
		syncIncidentFluxVec();

		return;
	}

	template <typename NetworkType>
	void
	addClusters(NetworkType& network)
	{
		// To scale the values to the default dose rate of 1.0e-6 dpa s-1
		double factor = getFluxAmplitude() * 1.0e6;
		auto omega = network.getAtomicVolume();

		// Set the flux index corresponding the the single helium cluster here
		typename NetworkType::Composition comp =
			NetworkType::Composition::zero();
		comp[NetworkType::Species::He] = 1;

		auto cluster = network.findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe single helium cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		incidentFluxVec.push_back(
			std::vector<double>(1, factor * hePPM * 1.0e-12 / omega));

		// Look for interstitial now
		comp[NetworkType::Species::He] = 0;
		comp[NetworkType::Species::I] = 1;
		cluster = network.findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe single interstitial cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		incidentFluxVec.push_back(std::vector<double>(1, factor * 1.49e-5));

		// Look for vacancies now
		comp[NetworkType::Species::I] = 0;
		comp[NetworkType::Species::V] = 1;
		cluster = network.findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe single vacancy cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		incidentFluxVec.push_back(std::vector<double>(1, factor * 9.91e-6));
		comp[NetworkType::Species::V] = 2;
		cluster = network.findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe double vacancy cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		incidentFluxVec.push_back(std::vector<double>(1, factor * 1.51e-6));
		comp[NetworkType::Species::V] = 3;
		cluster = network.findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe triple vacancy cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		incidentFluxVec.push_back(std::vector<double>(1, factor * 2.60e-7));
		comp[NetworkType::Species::V] = 4;
		cluster = network.findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe quadruple vacancy cluster is not present "
				"in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		incidentFluxVec.push_back(std::vector<double>(1, factor * 1.58e-7));
		comp[NetworkType::Species::V] = 5;
		cluster = network.findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nVacancy 5 cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		incidentFluxVec.push_back(std::vector<double>(1, factor * 6.29e-8));
		comp[NetworkType::Species::V] = 9;
		cluster = network.findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nVacancy 9 cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		incidentFluxVec.push_back(std::vector<double>(1, factor * 3.16e-8));

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
		auto ids = this->fluxIds;
		auto flux = this->incidentFlux;
		Kokkos::parallel_for(
			ids.size(), KOKKOS_LAMBDA(std::size_t i) {
				Kokkos::atomic_add(&updatedConcOffset[ids[i]], flux(i, 0));
			});
	}
};
// end class FeFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl
