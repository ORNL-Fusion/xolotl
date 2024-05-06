#ifndef ALLOYFLUXHANDLER_H
#define ALLOYFLUXHANDLER_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/AlloyReactionNetwork.h>
#include <xolotl/util/MPIUtils.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the FluxHandler interface to calculate the incident
 * fluxes for the alloy case.
 */
class AlloyFluxHandler : public FluxHandler
{
protected:
	/**
	 * \see FluxHandler.h
	 */
	double
	FitFunction(double x)
	{
		// Not actually used
		return 1.0;
	}

	std::vector<double> fluxI;
	std::vector<double> fluxV;

	double perfectFraction = 0.2;

public:
	/**
	 * The constructor
	 */
	AlloyFluxHandler(const options::IOptions& options) : FluxHandler(options)
	{
	}

	/**
	 * The Destructor
	 */
	~AlloyFluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid)
	{
		// Set the grid
		xGrid = grid;

		auto xolotlComm = util::getMPIComm();
		int procId;
		MPI_Comm_rank(xolotlComm, &procId);

		using NetworkType = network::AlloyReactionNetwork;
		auto alloyNetwork = dynamic_cast<NetworkType*>(&network);

		auto omega = alloyNetwork->getAtomicVolume();
		auto fluxFactor = fluxAmplitude / omega;

		// Set the flux index corresponding the interstitial clusters
		NetworkType::Composition comp = NetworkType::Composition::zero();
		for (int i = 0; i < fluxI.size(); i++) {
			comp[NetworkType::Species::I] = i;
			auto cluster =
				alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
			if (cluster.getId() == NetworkType::invalidIndex()) {
				continue;
			}
			fluxIndices.push_back(cluster.getId());
			std::vector<double> tempVector;
			if (xGrid.size() == 0)
				tempVector.push_back(fluxI[i] * fluxFactor);
			else {
				for (auto i = 0; i < xGrid.size(); i++) {
					tempVector.push_back(fluxI[i] * fluxFactor);
				}
			}
			incidentFluxVec.push_back(tempVector);
		}

		comp[NetworkType::Species::I] = 0;

		// Set the flux index corresponding the interstitial loops
		for (int i = 0; i < fluxI.size(); i++) {
			// Perfect
			comp[NetworkType::Species::FaultedI] = 0;
			comp[NetworkType::Species::PerfectI] = i;
			auto cluster =
				alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
			if (cluster.getId() == NetworkType::invalidIndex()) {
				continue;
			}
			fluxIndices.push_back(cluster.getId());
			std::vector<double> tempVector;
			if (xGrid.size() == 0)
				tempVector.push_back(fluxI[i] * perfectFraction * fluxFactor);
			else {
				for (auto i = 0; i < xGrid.size(); i++) {
					tempVector.push_back(
						fluxI[i] * perfectFraction * fluxFactor);
				}
			}
			incidentFluxVec.push_back(tempVector);

			// Faulted
			comp[NetworkType::Species::PerfectI] = 0;
			comp[NetworkType::Species::FaultedI] = i;
			cluster = alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
			if (cluster.getId() == NetworkType::invalidIndex()) {
				continue;
			}
			fluxIndices.push_back(cluster.getId());
			tempVector.clear();
			if (xGrid.size() == 0)
				tempVector.push_back(
					fluxI[i] * (1.0 - perfectFraction) * fluxFactor);
			else {
				for (auto i = 0; i < xGrid.size(); i++) {
					tempVector.push_back(
						fluxI[i] * (1.0 - perfectFraction) * fluxFactor);
				}
			}
			incidentFluxVec.push_back(tempVector);
		}

		comp[NetworkType::Species::FaultedI] = 0;
		comp[NetworkType::Species::PerfectI] = 0;

		// Set the flux index corresponding the vacancy clusters
		for (int i = 0; i < fluxV.size(); i++) {
			comp[NetworkType::Species::V] = i;
			auto cluster =
				alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
			if (cluster.getId() == NetworkType::invalidIndex()) {
				continue;
			}
			fluxIndices.push_back(cluster.getId());
			std::vector<double> tempVector;
			if (xGrid.size() == 0)
				tempVector.push_back(fluxV[i] * fluxFactor);
			else {
				for (auto i = 0; i < xGrid.size(); i++) {
					tempVector.push_back(fluxV[i] * fluxFactor);
				}
			}
			incidentFluxVec.push_back(tempVector);
		}

		comp[NetworkType::Species::V] = 0;

		// Set the flux index corresponding the vacancy loops
		for (int i = 0; i < fluxV.size(); i++) {
			// Perfect
			comp[NetworkType::Species::FaultedV] = 0;
			comp[NetworkType::Species::PerfectV] = i;
			auto cluster =
				alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
			if (cluster.getId() == NetworkType::invalidIndex()) {
				continue;
			}
			fluxIndices.push_back(cluster.getId());
			std::vector<double> tempVector;
			if (xGrid.size() == 0)
				tempVector.push_back(fluxV[i] * perfectFraction * fluxFactor);
			else {
				for (auto i = 0; i < xGrid.size(); i++) {
					tempVector.push_back(
						fluxV[i] * perfectFraction * fluxFactor);
				}
			}
			incidentFluxVec.push_back(tempVector);

			// Faulted
			comp[NetworkType::Species::PerfectV] = 0;
			comp[NetworkType::Species::FaultedV] = i;
			cluster = alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
			if (cluster.getId() == NetworkType::invalidIndex()) {
				continue;
			}
			fluxIndices.push_back(cluster.getId());
			tempVector.clear();
			if (xGrid.size() == 0)
				tempVector.push_back(
					fluxV[i] * (1.0 - perfectFraction) * fluxFactor);
			else {
				for (auto i = 0; i < xGrid.size(); i++) {
					tempVector.push_back(
						fluxV[i] * (1.0 - perfectFraction) * fluxFactor);
				}
			}
			incidentFluxVec.push_back(tempVector);
		}

		if (procId == 0) {
			std::ofstream outfile;
			outfile.open("alloyFlux.dat");
			for (int it = 0; it < fluxIndices.size(); ++it) {
				outfile << fluxIndices[it] << ": ";
				for (auto xi = 0; xi < std::max((int)grid.size(), 1); xi++) {
					outfile << incidentFluxVec[it][xi] << " ";
				}
				outfile << std::endl;
			}
			outfile.close();
		}

		// Sync data
		syncFluxIndices();
		syncIncidentFluxVec();
	}

	/**
	 * This operation computes the flux due to incoming particles at a given
	 * grid point. \see IFluxHandler.h
	 */
	void
	computeIncidentFlux(double currentTime,
		Kokkos::View<double*> updatedConcOffset, int xi,
		int surfacePos) override
	{
		// Attenuation factor to model reduced production of new point defects
		// with increasing dose (or time).
		double attenuation = 1.0;
		if (cascadeDose > 0.0) {
			attenuation = ((1.0 - cascadeEfficiency) / 2.0) *
					(1.0 -
						tanh(47.0 *
							(currentTime * fluxAmplitude - cascadeDose))) +
				cascadeEfficiency;
		}

		// Update the concentration array
		auto ids = this->fluxIds;
		auto flux = this->incidentFlux;
		Kokkos::parallel_for(
			ids.size(), KOKKOS_LAMBDA(std::size_t i) {
				Kokkos::atomic_add(
					&updatedConcOffset[ids[i]], attenuation * flux(i, xi));
			});
	}
}; // namespace flux
// end class AlloyFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
