#pragma once

#include <cmath>

#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/LiReactionNetwork.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the IFluxHandler interface to calculate the incident
 * fluxes for an lithium material.
 */
class LiFitFluxHandler : public FluxHandler
{
private:
	/**
	 * \see FluxHandler.h
	 */
	double
	FitFunction(double x) override
	{
		// Flat profile
		return 1.0;
	}

public:
	/**
	 * The constructor
	 */
	LiFitFluxHandler(const options::IOptions& options) : FluxHandler(options)
	{
	}

	/**
	 * The Destructor
	 */
	~LiFitFluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid) override
	{
		// Call the general method
		FluxHandler::initializeFluxHandler(network, surfacePos, grid);

		// Skip if the flux amplitude is 0.0 and we are not using a time profile
		if (util::equal(fluxAmplitude, 0.0) && !useTimeProfile)
			return;

		using NetworkType = network::LiReactionNetwork;
		auto liNetwork = dynamic_cast<NetworkType*>(&network);

		// Set the flux index corresponding the the single helium cluster here
		NetworkType::Composition comp = NetworkType::Composition::zero();
		comp[NetworkType::Species::H] = 1;
		auto cluster = liNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe single hydrogen cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		
		// Vacancy
		comp[NetworkType::Species::H] = 0;
		comp[NetworkType::Species::V] = 1;
		cluster = liNetwork->findCluster(comp, plsm::HostMemSpace{});
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe single vacancy cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		return;
	}

	/**
	 * \see IFluxHandler.h
	 */
        void
        computeIncidentFlux(double currentTime,
	        Kokkos::View<const double*>, Kokkos::View<double*> updatedConcOffset,
	        int xi, int surfacePos)
        {
	        // Skip if no index was set
	        if (fluxIndices.size() == 0)
		        return;

	        // Recompute the flux vector if a time profile is used
	        if (useTimeProfile) {
		        fluxAmplitude = getProfileAmplitude(currentTime);
		        recomputeFluxHandler(surfacePos);
	        }

	        double value{};
	        if (incidentFluxVec[0].size() == 0) {
		        value = fluxAmplitude;
	        }
	        else {
		        value = incidentFluxVec[0][xi - surfacePos];
	        }

		auto ids = fluxIndices;
	        // Update the concentration array
	        Kokkos::parallel_for(
		        ids.size(), KOKKOS_LAMBDA(std::size_t i) { updatedConcOffset[ids[i]] += value; });
        }
};
// end class LiFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl
