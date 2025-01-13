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
	std::vector<double> highFluxI;
	std::vector<double> fluxV;

	double perfectFraction = 0.2;
	double voidFraction = 0.5;

	/**
	 * Vector to hold the incident flux values at each grid
	 * point (x position).
	 */
	std::vector<std::vector<double>> incidentHighFluxVec;
	Kokkos::View<double**> incidentHighFlux;

	/**
	 * The indices of the incoming clusters.
	 */
	std::vector<IdType> highFluxIndices;
	Kokkos::View<IdType*> highFluxIds;

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

		// Set the flux index corresponding the high energy ions interstitial
		// clusters
		for (int i = 0; i < highFluxI.size(); i++) {
			comp[NetworkType::Species::I] = i;
			auto cluster =
				alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
			if (cluster.getId() == NetworkType::invalidIndex()) {
				continue;
			}
			fluxIndices.push_back(cluster.getId());
			std::vector<double> tempVector;
			if (xGrid.size() == 0)
				tempVector.push_back(highFluxI[i] * fluxFactor);
			else {
				for (auto i = 0; i < xGrid.size(); i++) {
					tempVector.push_back(highFluxI[i] * fluxFactor);
				}
			}
			incidentFluxVec.push_back(tempVector);
		}

		comp[NetworkType::Species::I] = 0;

		// Set the flux index corresponding the high energy ions interstitial
		// loops
		for (int i = 0; i < highFluxI.size(); i++) {
			// Perfect
			comp[NetworkType::Species::FaultedI] = 0;
			comp[NetworkType::Species::PerfectI] = i;
			auto cluster =
				alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
			if (cluster.getId() == NetworkType::invalidIndex()) {
				continue;
			}
			highFluxIndices.push_back(cluster.getId());
			std::vector<double> tempVector;
			if (xGrid.size() == 0)
				tempVector.push_back(
					highFluxI[i] * perfectFraction * fluxFactor);
			else {
				for (auto i = 0; i < xGrid.size(); i++) {
					tempVector.push_back(
						highFluxI[i] * perfectFraction * fluxFactor);
				}
			}
			incidentHighFluxVec.push_back(tempVector);

			// Faulted
			comp[NetworkType::Species::PerfectI] = 0;
			comp[NetworkType::Species::FaultedI] = i;
			cluster = alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
			if (cluster.getId() == NetworkType::invalidIndex()) {
				continue;
			}
			highFluxIndices.push_back(cluster.getId());
			tempVector.clear();
			if (xGrid.size() == 0)
				tempVector.push_back(
					highFluxI[i] * (1.0 - perfectFraction) * fluxFactor);
			else {
				for (auto i = 0; i < xGrid.size(); i++) {
					tempVector.push_back(
						highFluxI[i] * (1.0 - perfectFraction) * fluxFactor);
				}
			}
			incidentHighFluxVec.push_back(tempVector);
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

			// Check if it is mobile or not for voids
			if (cluster.getDiffusionFactor() == 0.0) {
				// Void
				if (xGrid.size() == 0)
					tempVector.push_back(fluxV[i] * voidFraction * fluxFactor);
				else {
					for (auto i = 0; i < xGrid.size(); i++) {
						tempVector.push_back(
							fluxV[i] * voidFraction * fluxFactor);
					}
				}
			}
			else {
				// Vacancy
				if (xGrid.size() == 0)
					tempVector.push_back(fluxV[i] * fluxFactor);
				else {
					for (auto i = 0; i < xGrid.size(); i++) {
						tempVector.push_back(fluxV[i] * fluxFactor);
					}
				}
			}
			incidentFluxVec.push_back(tempVector);
		}

		comp[NetworkType::Species::V] = 0;

		// Set the flux index corresponding the vacancy loops
		for (int i = 0; i < fluxV.size(); i++) {
			// Faulted
			comp[NetworkType::Species::FaultedV] = i;
			auto cluster =
				alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
			if (cluster.getId() == NetworkType::invalidIndex()) {
				continue;
			}
			fluxIndices.push_back(cluster.getId());
			std::vector<double> tempVector;
			if (xGrid.size() == 0)
				tempVector.push_back(
					fluxV[i] * (1.0 - voidFraction) * fluxFactor);
			else {
				for (auto i = 0; i < xGrid.size(); i++) {
					tempVector.push_back(
						fluxV[i] * (1.0 - voidFraction) * fluxFactor);
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

		// Sync high energy data
		auto ids_h =
			Kokkos::View<IdType*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
				highFluxIndices.data(), highFluxIndices.size());
		highFluxIds = Kokkos::View<IdType*>(
			Kokkos::ViewAllocateWithoutInitializing("Flux Indices"),
			highFluxIndices.size());
		deep_copy(highFluxIds, ids_h);

		incidentHighFlux = Kokkos::View<double**>("Incident High Flux Vec",
			incidentHighFluxVec.size(), incidentHighFluxVec[0].size());
		auto incidentFlux_h = create_mirror_view(incidentHighFlux);
		for (std::size_t i = 0; i < incidentHighFluxVec.size(); ++i) {
			for (std::size_t j = 0; j < incidentHighFluxVec[i].size(); ++j) {
				incidentFlux_h(i, j) = incidentHighFluxVec[i][j];
			}
		}
		deep_copy(incidentHighFlux, incidentFlux_h);
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
		// Update the concentration array
		auto ids = this->fluxIds;
		auto flux = this->incidentFlux;
		Kokkos::parallel_for(
			ids.size(), KOKKOS_LAMBDA(std::size_t i) {
				Kokkos::atomic_add(&updatedConcOffset[ids[i]], flux(i, xi));
			});
		// Update the concentration array for high energy ions
		ids = this->highFluxIds;
		flux = this->incidentHighFlux;
		Kokkos::parallel_for(
			ids.size(), KOKKOS_LAMBDA(std::size_t i) {
				Kokkos::atomic_add(&updatedConcOffset[ids[i]],
					(1.0 - deltaCorrection) * flux(i, xi));
			});
	}

	/**
	 * \see IFluxHandler.h
	 */
	std::vector<double>
	getHighFluxVector() const
	{
		return highFluxI;
	}
}; // namespace flux
// end class AlloyFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
