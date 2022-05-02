#ifndef FECRFITFLUXHANDLER_H
#define FECRFITFLUXHANDLER_H

#include <cmath>

#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/FeCrReactionNetwork.h>

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
class FeCrFitFluxHandler : public FluxHandler
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
	FeCrFitFluxHandler(const options::IOptions& options) : FluxHandler(options)
	{
	}

	/**
	 * The Destructor
	 */
	~FeCrFitFluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid)
	{
		using NetworkType = network::FeCrReactionNetwork;
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

		// Free
		comp[NetworkType::Species::I] = 0;
		comp[NetworkType::Species::Free] = 4;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe Free4 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::Free] = 5;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe Free5 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::Free] = 6;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe Free6 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::Free] = 7;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe Free7 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::Free] = 8;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe Free8 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::Free] = 9;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe Free9 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::Free] = 12;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe Free12 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::Free] = 16;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe Free16 cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::Free] = 20;
		cluster = feNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe Free20 cluster is not "
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
	computeIncidentFlux(
		double currentTime, double* updatedConcOffset, int xi, int surfacePos)
	{
		updatedConcOffset[fluxIndices[0]] += fluxAmplitude * 11.7769; // V1
		updatedConcOffset[fluxIndices[1]] += fluxAmplitude * 1.815; // V2
		updatedConcOffset[fluxIndices[2]] += fluxAmplitude * 0.3146; // V3
		updatedConcOffset[fluxIndices[3]] += fluxAmplitude * 0.1914; // V4
		updatedConcOffset[fluxIndices[4]] += fluxAmplitude * 0.07611; // V5
		updatedConcOffset[fluxIndices[5]] += fluxAmplitude * 0.03826; // V9
		updatedConcOffset[fluxIndices[6]] += fluxAmplitude * 5.94; // I1
		updatedConcOffset[fluxIndices[7]] += fluxAmplitude * 2.1869; // I2
		updatedConcOffset[fluxIndices[8]] += fluxAmplitude * 0.7003; // I3
		updatedConcOffset[fluxIndices[9]] += fluxAmplitude * 0.4478; // Free4
		updatedConcOffset[fluxIndices[10]] += fluxAmplitude * 0.2619; // Free5
		updatedConcOffset[fluxIndices[11]] += fluxAmplitude * 0.1615; // Free6
		updatedConcOffset[fluxIndices[12]] += fluxAmplitude * 0.09497; // Free7
		updatedConcOffset[fluxIndices[13]] += fluxAmplitude * 0.02311; // Free8
		updatedConcOffset[fluxIndices[14]] += fluxAmplitude * 0.02581; // Free9
		updatedConcOffset[fluxIndices[15]] +=
			fluxAmplitude * 0.009665; // Free12
		updatedConcOffset[fluxIndices[16]] +=
			fluxAmplitude * 0.006116; // Free16
		updatedConcOffset[fluxIndices[17]] +=
			fluxAmplitude * 0.003081; // Free20

		return;
	}
};
// end class FeCrFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
