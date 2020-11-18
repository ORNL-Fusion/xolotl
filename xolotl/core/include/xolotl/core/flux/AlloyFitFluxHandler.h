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
	//! The time parameter for attenuation
	double tauFlux;

	/**
	 * \see FluxHandler.h
	 */
	double
	FitFunction(double x)
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

	struct IonDamageStruct
	{
		std::vector<int> fluxIndex;
		std::vector<std::vector<double>> damageRate;
	} ionDamage;

	Cascade cascade;
	SRIMData srim;

public:
	/**
	 * The constructor
	 */
	AlloyFitFluxHandler(const options::IOptions& options) :
		FluxHandler(options),
		tauFlux(0.0)
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
					alloyNetwork->findCluster(comp, plsm::onHost);
				if (fluxCluster.getId() != NetworkType::invalidIndex()) {
					(ionDamage.fluxIndex).push_back(fluxCluster.getId());
					(ionDamage.damageRate)
						.push_back(AlloySetGeneration(size, it, 1.0));
					if (size == 1 && implant) {
						AlloyAddImplantation(ionDamage.damageRate.back());
					}
				}
				// Otherwise the clusters must be frank and perfect type
				else {
					comp[NetworkType::Species::I] = 0;
					comp[NetworkType::Species::Frank] = size;
					auto fluxCluster1 =
						alloyNetwork->findCluster(comp, plsm::onHost);
					comp[NetworkType::Species::Frank] = 0;
					comp[NetworkType::Species::Perfect] = size;
					auto fluxCluster2 =
						alloyNetwork->findCluster(comp, plsm::onHost);
					if (fluxCluster1.getId() == NetworkType::invalidIndex() ||
						fluxCluster2.getId() == NetworkType::invalidIndex()) {
						// Throw error -> missing type
						throw std::runtime_error(
							"\nNo clusted of size: " + std::to_string(size) +
							", cannot use the flux option!");
					}
					else {
						// Frank loop
						(ionDamage.fluxIndex).push_back(fluxCluster1.getId());
						double frac = 1.0 - cascade.perfectFraction;
						(ionDamage.damageRate)
							.push_back(AlloySetGeneration(size, it, frac));
						// Perfect loop
						(ionDamage.fluxIndex).push_back(fluxCluster2.getId());
						frac = cascade.perfectFraction;
						(ionDamage.damageRate)
							.push_back(AlloySetGeneration(size, it, frac));
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
					alloyNetwork->findCluster(comp, plsm::onHost);
				if (fluxCluster.getId() != NetworkType::invalidIndex()) {
					(ionDamage.fluxIndex).push_back(fluxCluster.getId());
					(ionDamage.damageRate)
						.push_back(AlloySetGeneration(size, it, 1.0));
				}
				// Otherwise the clusters must be faulted type
				else {
					comp[NetworkType::Species::V] = 0;
					comp[NetworkType::Species::Faulted] = size;
					fluxCluster = alloyNetwork->findCluster(comp, plsm::onHost);
					if (fluxCluster.getId() == NetworkType::invalidIndex()) {
						// Throw error -> no available type
						throw std::runtime_error(
							"\nNo clusted of size: " + std::to_string(-size) +
							", cannot use the flux option!");
					}
					else {
						// Faulted loop
						(ionDamage.fluxIndex).push_back(fluxCluster.getId());
						(ionDamage.damageRate)
							.push_back(AlloySetGeneration(size, it, 1.0));
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

		if (procId == 0) {
			std::ofstream outfile;
			outfile.open("alloyFlux.dat");
			for (int it = 0; it < ionDamage.fluxIndex.size(); ++it) {
				outfile << ionDamage.fluxIndex[it] << ": ";
				for (int xi = surfacePos; xi < std::max((int)grid.size(), 1);
					 xi++) {
					outfile << ionDamage.damageRate[it][xi - surfacePos] << " ";
				}
				outfile << std::endl;
			}
			outfile.close();
		}

		return;
	}

	/**
	 * This operation computes the flux due to incoming particles at a given
	 * grid point. \see IFluxHandler.h
	 */
	void
	computeIncidentFlux(
		double currentTime, double* updatedConcOffset, int xi, int surfacePos)
	{
		// Attenuation factor to model reduced production of new point defects
		// with increasing dose (or time).
		double attenuation = 1.0;
		if (tauFlux > 0.0 && currentTime > 0.0)
			attenuation = 1.0 - exp((-1.0 * tauFlux) / currentTime);

		// Update the concentration array
		for (int it = 0; it < ionDamage.fluxIndex.size(); ++it) {
			updatedConcOffset[ionDamage.fluxIndex[it]] +=
				attenuation * ionDamage.damageRate[it][xi - surfacePos];
		}

		return;
	}

	/**
	 * This operation sets the attenuation parameter.
	 * \see IFluxHandler.h
	 */
	void
	setTauFlux(double tau)
	{
		tauFlux = tau;
		return;
	}
};
// end class AlloyFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
