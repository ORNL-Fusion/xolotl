#ifndef ALLOYFITFLUXHANDLER_H
#define ALLOYFITFLUXHANDLER_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include <xolotl/core/flux/AlloySRIMData.h>
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
class AlloyFitFluxHandler : public FluxHandler
{
private:
	/**
	 * \see FluxHandler.h
	 */
	double
	FitFunction(double x)
	{
		// Not actually used
		return 1.0;
	}

	std::vector<double> fluxI = {0.0, 618.4963326, 3.476762215, 1.408162213,
		0.982739683, 0.30225657, 0.327148562, 0.3094983, 0.12474936,
		0.088058372, 0.133219224, 0.0, 0.0, 0.09539657, 0.0, 0.0, 0.0, 0.0,
		0.102734767, 0.0, 0.0, 0.0, 0.0, 0.09539657, 0.0, 0.0, 0.0, 0.0,
		0.058705581, 0.0, 0.0, 0.0, 0.0, 0.011007296, 0.0, 0.0, 0.0, 0.0,
		0.011007297, 0.0, 0.0, 0.0, 0.0, 0.011007297};

	std::vector<double> fluxV = {0.0, 634.8181526, 1.05768736, 0.320601267,
		0.249215543, 0.16824146, 0.088914908, 0.17132574, 0.050966698,
		0.054007957, 0.052681762, 0.0, 0.0, 0.094279052, 0.0, 0.0, 0.0, 0.0,
		0.022809448, 0.0, 0.0, 0.0, 0.0, 0.045618896, 0.0, 0.0, 0.0, 0.0,
		0.019768188, 0.0, 0.0, 0.0, 0.0, 0.027371338, 0.0, 0.0, 0.0, 0.0,
		0.013685669, 0.0, 0.0, 0.0, 0.0, 0.009123779, 0.0, 0.0, 0.0, 0.0,
		0.021288818, 0.0, 0.0, 0.0, 0.0, 0.012165039, 0.0, 0.0, 0.0, 0.0,
		0.009123779, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00608252};

	double perfectFraction = 0.2;

public:
	/**
	 * The constructor
	 */
	AlloyFitFluxHandler(const options::IOptions& options) : FluxHandler(options)
	{
	}

	/**
	 * The Destructor
	 */
	~AlloyFitFluxHandler()
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
// end class AlloyFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
