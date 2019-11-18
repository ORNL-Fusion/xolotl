#ifndef ALLOYFITFLUXHANDLER_H
#define ALLOYFITFLUXHANDLER_H

#include "FluxHandler.h"
#include "AlloySRIMData.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident fluxes
 * for the alloy case.
 */
class AlloyFitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 *
	 * @param x The position where to evaluate he fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		return 1.0;
	}

	std::vector<double> AlloyDamageFunction(const double x) {
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

	double AlloyImplantationFunction(const double x) {
		// Find the correct depth region
		for (int i = 0; i < srim.getDepth().size(); ++i) {
			if (x <= srim.getDepth()[i]) {
				return srim.getImplantation()[i];
			}
		}
		return 0.0;
	}

	std::vector<double> AlloySetGeneration(const int size, const int it,
			const double fraction) {
		std::vector<double> damageRate;
		// 0D case
		if (xGrid.size() == 0) {
			// Compute the rate at that position
			std::vector<double> fitFlux = AlloyDamageFunction(0);
			double rate = 0;
			for (int j = 0; j < fitFlux.size(); ++j) {
				rate += cascade.cascadeEfficiency[j] * fitFlux[j]
						* cascade.clusterFraction[it][j];
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
				rate += cascade.cascadeEfficiency[j] * fitFlux[j]
						* cascade.clusterFraction[it][j];
			}
			rate = rate * fraction * fluxAmplitude / double(size);

			// Add it to the vector
			damageRate.push_back(rate);
		}
		damageRate.push_back(0.0);
		return damageRate;
	}

	void AlloyAddImplantation(std::vector<double> & input) {
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

	struct IonDamageStruct {
		std::vector<int> fluxIndex;
		std::vector<std::vector<double> > damageRate;
	} ionDamage;

	Cascade cascade;
	SRIMData srim;

public:

	/**
	 * The constructor
	 */
	AlloyFitFluxHandler() {
	}

	/**
	 * The Destructor
	 */
	~AlloyFitFluxHandler() {
	}

	/**
	 * Compute and store the incident flux values at each grid point.
	 * \see IFluxHandler.h
	 */
	void initializeFluxHandler(const IReactionNetwork& network, int surfacePos,
			std::vector<double> grid) {

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

		// Iterate over all produced cluster species
		for (int it = 0; it < cascade.clusterSizes.size(); ++it) {

			// Get the size of the cluster
			int size = cascade.clusterSizes[it];

			// Check if cluster is interstitial type
			if (size > 0) {
				// See if theres an iType cluster of size
				auto fluxCluster = network.get(toSpecies(ReactantType::I),
						size);
				if (fluxCluster) {
					(ionDamage.fluxIndex).push_back(fluxCluster->getId() - 1);
					(ionDamage.damageRate).push_back(
							AlloySetGeneration(size, it, 1.0));
					if (size == 1 && implant) {
						AlloyAddImplantation(ionDamage.damageRate.back());
					}
				}
				// Otherwise the clusters must be frank and perfect type
				else {
					auto fluxCluster1 = network.get(
							toSpecies(ReactantType::Frank), size);
					auto fluxCluster2 = network.get(
							toSpecies(ReactantType::Perfect), size);
					if (!fluxCluster1 || !fluxCluster2) {
						// Throw error -> missing type
						std::cout << "Error: no flux cluster of size " << size
								<< std::endl;
					} else {
						// Frank loop
						(ionDamage.fluxIndex).push_back(
								fluxCluster1->getId() - 1);
						double frac = 1.0 - cascade.perfectFraction;
						(ionDamage.damageRate).push_back(
								AlloySetGeneration(size, it, frac));
						// Perfect loop
						(ionDamage.fluxIndex).push_back(
								fluxCluster2->getId() - 1);
						frac = cascade.perfectFraction;
						(ionDamage.damageRate).push_back(
								AlloySetGeneration(size, it, frac));
					}
				}
			}
			// Check if cluster is vacancy type
			else if (size < 0) {
				size = -size;
				// See if theres an vType cluster of size
				auto fluxCluster = network.get(toSpecies(ReactantType::V),
						size);
				if (fluxCluster) {
					(ionDamage.fluxIndex).push_back(fluxCluster->getId() - 1);
					(ionDamage.damageRate).push_back(
							AlloySetGeneration(size, it, 1.0));
				}
				// Otherwise the clusters must be faulted type
				else {
					fluxCluster = network.get(toSpecies(ReactantType::Faulted),
							size);
					if (!fluxCluster) {
						// Throw error -> no available type
						std::cout << "Error: no flux cluster of size " << -size
								<< std::endl;
					} else {
						// Faulted loop
						(ionDamage.fluxIndex).push_back(
								fluxCluster->getId() - 1);
						(ionDamage.damageRate).push_back(
								AlloySetGeneration(size, it, 1.0));
					}
				}
			}
			// Neither interstitial nor vacancy type cluster
			else {
				// Throw error for size 0 cluster
			}

		}

		// record the results in a files
		std::ofstream outfile;
		outfile.open("alloyFlux.dat");
		auto clusters = network.getAll();
		for (int i = 0; i < ionDamage.fluxIndex.size(); ++i) {
			outfile << ionDamage.fluxIndex[i];
			for (int j = 0; j < ionDamage.damageRate[i].size(); ++j) {
				outfile << " " << ionDamage.damageRate[i][j];
			}
			// find the corresponding flux clusters
			std::for_each(clusters.begin(), clusters.end(),
					[&outfile, &i, this](IReactant& currReactant) {
						if (currReactant.getId() - 1 == ionDamage.fluxIndex[i]) {
							outfile << " " << currReactant.getName();
						}
					});
			outfile << std::endl;
		}
		outfile.close();

		return;
	}

	/**
	 * This operation computes the flux due to incoming particles at a given grid point.
	 * \see IFluxHandler.h
	 */
	void computeIncidentFlux(double currentTime, double *updatedConcOffset,
			int xi, int surfacePos) {
		
		// Attenuation factor to model reduced production of new point defects
		// with increasing dose (or time).
		double tau = 1.0e7;
		double attenuation = 1.0 - exp((-1.0 * tau) / currentTime);

		// Update the concentration array
		for (int it = 0; it < ionDamage.fluxIndex.size(); ++it) {
			updatedConcOffset[ionDamage.fluxIndex[it]] +=
					attenuation * ionDamage.damageRate[it][xi - surfacePos];
		}

		return;
	}

};
// end class AlloyFitFluxHandler

}

#endif
