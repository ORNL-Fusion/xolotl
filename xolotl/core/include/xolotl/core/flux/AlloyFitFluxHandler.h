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
	FitFunction(double x) override
	{
		// Not actually used
		return 1.0;
	}

	/**
	 * Computes the damage as a function of the data from SRIM
	 *
	 * @param x The position
	 * @return The vector of damages
	 */
	std::vector<double>
	AlloyDamageFunction(const double x)
	{
		std::vector<double> damage;
		auto srimDamage = srim.getDamage();
		for (int it = 0; it < srimDamage.size(); ++it) {
			damage.push_back(srimDamage[it][0]);
			for (int it2 = 1; it2 < srimDamage[it].size(); ++it2) {
				damage[it] += srimDamage[it][it2] * pow(x, double(it2));
			}
		}
		return damage;
	}

	/**
	 * Computes the implantation as a function of the data from SRIM
	 *
	 * @param x The position
	 * @return The value of implantation
	 */
	double
	AlloyImplantationFunction(const double x)
	{
		// Find the correct depth region
		for (int i = 0; i < srim.getDepth().size(); ++i) {
			if (x <= srim.getDepth()[i]) {
				return srim.getImplantation()[i];
			}
		}
		return 0.0;
	}

	/**
	 * Computes the generation rate as a function of the data from SRIM
	 *
	 * @param size The size of the cluster
	 * @param it The index of the cluster in the cascade
	 * @param fraction The reduction fraction
	 * @return The vector of damage rates
	 */
	std::vector<double>
	AlloySetGeneration(const int size, const int it, const double fraction)
	{
		std::vector<double> damageRate;
		// 0D case
		if (xGrid.size() == 0) {
			// Compute the rate at that position
			std::vector<double> fitFlux = AlloyDamageFunction(0);
			double rate = 0;
			for (int j = 0; j < fitFlux.size(); ++j) {
				rate += cascade.cascadeEfficiency[j] * fitFlux[j] *
					cascade.clusterFraction[it][j];
			}
			rate = rate * fraction * fluxAmplitude / double(size);
			// Add it to the vector
			damageRate.push_back(rate);
			return damageRate;
		}

		// 1D case
		// Change this to grab the actual surface position
		int surfacePos = 0;
		damageRate.push_back(0.0);
		for (int i = surfacePos + 1; i < xGrid.size() - 3; i++) {
			// Get the x position
			auto x = xGrid[i + 1] - xGrid[surfacePos + 1];
			// Compute the rate at that position
			std::vector<double> fitFlux = AlloyDamageFunction(x);
			double rate = 0;
			for (int j = 0; j < fitFlux.size(); ++j) {
				rate += cascade.cascadeEfficiency[j] * fitFlux[j] *
					cascade.clusterFraction[it][j];
			}
			rate = rate * fraction * fluxAmplitude / double(size);

			// Add it to the vector
			damageRate.push_back(rate);
		}
		damageRate.push_back(0.0);
		return damageRate;
	}

	/**
	 * Adds contribution to the given input
	 *
	 * @param input The input to add to
	 */
	void
	AlloyAddImplantation(std::vector<double>& input)
	{
		// Change this to grab the actual surface position
		int surfacePos = 0;

		for (int i = surfacePos + 1; i < xGrid.size() - 3; i++) {
			// Get the x position
			auto x = xGrid[i] - xGrid[surfacePos];
			// Add the implantation rate to the damage rate
			input[i] += AlloyImplantationFunction(x);
		}

		return;
	}

	struct IonDamage
	{
		Kokkos::View<IdType*> fluxIds;
		Kokkos::View<double**> rate;
	} ionDamage;

	Cascade cascade;
	SRIMData srim;

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
		std::vector<double> grid) override
	{
		// Setup the ion damage and implantation depth profile
		if (false) {
			srim.setInSitu();
			cascade.setBulk();
		}
		else if (true) {
			srim.setBulk();
			cascade.setBulk();
		}

		// Turn on/off implantation
		bool implant = false;

		// Set the grid
		xGrid = grid;

		if (xGrid.size() == 0) {
			srim.setOverlap();
			cascade.setOverlap();
		}
		auto xolotlComm = util::getMPIComm();
		int procId;
		MPI_Comm_rank(xolotlComm, &procId);

		using NetworkType = network::AlloyReactionNetwork;
		auto alloyNetwork = dynamic_cast<NetworkType*>(&network);

		std::vector<IdType> damageIds;
		std::vector<std::vector<double>> damageRates;

		// Iterate over all produced cluster species
		for (int it = 0; it < cascade.clusterSizes.size(); ++it) {
			// Get the size of the cluster
			int size = cascade.clusterSizes[it];

			// Check if cluster is interstitial type
			if (size > 0) {
				// See if theres an iType cluster of size
				NetworkType::Composition comp =
					NetworkType::Composition::zero();
				comp[NetworkType::Species::I] = size;
				auto fluxCluster =
					alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
				if (fluxCluster.getId() != NetworkType::invalidIndex()) {
					damageIds.push_back(fluxCluster.getId());
					damageRates.push_back(AlloySetGeneration(size, it, 1.0));
					if (size == 1 && implant) {
						AlloyAddImplantation(damageRates.back());
					}
				}
				// Otherwise the clusters must be frank and perfect type
				else {
					comp[NetworkType::Species::I] = 0;
					comp[NetworkType::Species::Frank] = size;
					auto fluxCluster1 =
						alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
					comp[NetworkType::Species::Frank] = 0;
					comp[NetworkType::Species::Perfect] = size;
					auto fluxCluster2 =
						alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
					if (fluxCluster1.getId() == NetworkType::invalidIndex() ||
						fluxCluster2.getId() == NetworkType::invalidIndex()) {
						continue;
					}
					else {
						// Frank loop
						damageIds.push_back(fluxCluster1.getId());
						double frac = 1.0 - cascade.perfectFraction;
						damageRates.push_back(
							AlloySetGeneration(size, it, frac));
						// Perfect loop
						damageIds.push_back(fluxCluster2.getId());
						frac = cascade.perfectFraction;
						damageRates.push_back(
							AlloySetGeneration(size, it, frac));
					}
				}
			}
			// Check if cluster is vacancy type
			else if (size < 0) {
				size = -size;
				// See if theres an vType cluster of size
				NetworkType::Composition comp =
					NetworkType::Composition::zero();
				comp[NetworkType::Species::V] = size;
				auto fluxCluster =
					alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
				if (fluxCluster.getId() != NetworkType::invalidIndex()) {
					damageIds.push_back(fluxCluster.getId());
					damageRates.push_back(AlloySetGeneration(size, it, 1.0));
				}
				// Otherwise the clusters must be faulted type
				else {
					comp[NetworkType::Species::V] = 0;
					comp[NetworkType::Species::Faulted] = size;
					fluxCluster =
						alloyNetwork->findCluster(comp, plsm::HostMemSpace{});
					if (fluxCluster.getId() == NetworkType::invalidIndex()) {
						continue;
					}
					else {
						// Faulted loop
						damageIds.push_back(fluxCluster.getId());
						damageRates.push_back(
							AlloySetGeneration(size, it, 1.0));
					}
				}
			}
			// Neither interstitial nor vacancy type cluster
			else {
				// Throw error for size 0 cluster
				throw std::runtime_error(
					"\nThe cluster is of size 0 which is not possible, cannot "
					"use the flux option!");
			}
		}

		if (damageIds.size() != damageRates.size()) {
			throw std::runtime_error("Ion damage ids and rates should have the "
									 "same number of entries.");
		}

		if (procId == 0) {
			std::ofstream outfile;
			outfile.open("alloyFlux.dat");
			for (int it = 0; it < damageIds.size(); ++it) {
				outfile << damageIds[it] << ": ";
				for (int xi = surfacePos; xi < std::max((int)grid.size(), 1);
					 xi++) {
					outfile << damageRates[it][xi - surfacePos] << " ";
				}
				outfile << std::endl;
			}
			outfile.close();
		}

		// Move ion damage data to device
		std::size_t nDamageVals = damageIds.size();
		Kokkos::View<IdType*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
			ionDamageFluxIds_h(damageIds.data(), nDamageVals);
		auto innerSize = xGrid.size() == 0 ? 1 : xGrid.size() - surfacePos + 1;
		ionDamage = IonDamage{
			Kokkos::View<IdType*>{"Ion Damage Flux Indices", nDamageVals},
			Kokkos::View<double**>{"Ion Damage Rate", nDamageVals, innerSize}};
		deep_copy(ionDamage.fluxIds, ionDamageFluxIds_h);

		auto ionDamageRate_h = create_mirror_view(ionDamage.rate);
		for (std::size_t i = 0; i < nDamageVals; ++i) {
			for (std::size_t j = 0; j < damageRates[i].size(); ++j) {
				ionDamageRate_h(i, j) = damageRates[i][j];
			}
		}
		deep_copy(ionDamage.rate, ionDamageRate_h);
	}

	/**
	 * This operation computes the flux due to incoming particles at a given
	 * grid point. \see IFluxHandler.h
	 */
	void
	computeIncidentFlux(double currentTime, Kokkos::View<const double*>,
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
		auto ionDamageFluxIds = this->ionDamage.fluxIds;
		auto ionDamageRate = this->ionDamage.rate;
		Kokkos::parallel_for(
			ionDamageFluxIds.size(), KOKKOS_LAMBDA(std::size_t i) {
				Kokkos::atomic_add(&updatedConcOffset[ionDamageFluxIds[i]],
					attenuation * ionDamageRate(i, xi - surfacePos));
			});
	}

	/**
	 * \see IFluxHandler.h
	 */
	std::vector<std::pair<IdType, double>>
	getImplantedFlux(std::vector<IdType> map) override
	{
		std::vector<std::pair<IdType, double>> toReturn;

		auto ionDamageFluxIds_h = create_mirror_view(ionDamage.fluxIds);
		auto ionDamageRate_h = create_mirror_view(ionDamage.rate);
		deep_copy(ionDamageFluxIds_h, ionDamage.fluxIds);
		deep_copy(ionDamageRate_h, ionDamage.rate);
		for (auto i = 0; i < map.size(); i++) {
			// Look for this value in fluxIndices
			for (auto j = 0; j < ionDamageFluxIds_h.size(); j++) {
				if (map[i] == ionDamageFluxIds_h[j]) {
					toReturn.push_back(
						std::make_pair(i, ionDamageRate_h(j, 0)));
					break;
				}
			}
		}
		return toReturn;
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	setImplantedFlux(std::vector<std::pair<IdType, double>> fluxVector) override
	{
		std::size_t nDamageVals = ionDamage.fluxIds.size();
		if (fluxVector.size() != nDamageVals) {
			throw std::runtime_error(
				"AlloyFitFluxHandler::setImplantedFlux: "
				"called with different number of ion damage values than "
				"determined during initialization.");
		}

		// Loop on the flux vector
		auto ionDamageFluxIds_h = create_mirror_view(ionDamage.fluxIds);
		auto ionDamageRate_h = create_mirror_view(ionDamage.rate);
		for (auto i = 0; i < nDamageVals; i++) {
			ionDamageFluxIds_h(i) = fluxVector[i].first;
			ionDamageRate_h(i, 0) = fluxVector[i].second;
		}

		deep_copy(ionDamage.fluxIds, ionDamageFluxIds_h);
		deep_copy(ionDamage.rate, ionDamageRate_h);
	}
};
// end class AlloyFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
