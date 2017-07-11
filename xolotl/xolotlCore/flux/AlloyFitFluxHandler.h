#ifndef ALLOYFITFLUXHANDLER_H
#define ALLOYFITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>
#include <iostream>

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

  /**
   * Function that calculates the displacement rate in disp/ion/nm at a given
	 * position x in nm and for different PKA energy ranges.
	 */
	std::vector<double> AlloyFitFunction (const double x) {
		std::vector<double> damage = {0.0,0.0,0.0,0.0};
		if (x > ionDamage.fitRange[0] && x < ionDamage.fitRange[1]) {
			for (int it = 0; it < ionDamage.fitCoefficients.size(); ++it) {
				damage[it] = ionDamage.fitCoefficients[it][0];
				for (int it2 = 1; it2 < ionDamage.fitCoefficients[it].size(); ++it2) {
					damage[it] += ionDamage.fitCoefficients[it][it2] *
							pow(x,double(it2));
				}
			}
		}
		return damage;
	}

  /**
   * Function that calculates the displacement rate for given cluster size (size)
	 * and identification (it) for all positions.
	 */
	std::vector<double> AlloySetDamage (const int size, const int it,
			const double fraction) {

		// Not sure how to grab surfacePos from other files
		int surfacePos = 0;
		std::vector<double> damageRate;
		damageRate.push_back(0.0);
		for (int i = surfacePos + 1; i < xGrid.size() - 1; i++) {
			// Get the x position
			auto x = xGrid[i] - xGrid[surfacePos];
			// Compute the rate at that position
			std::vector<double> fitFlux = AlloyFitFunction(x);
			double rate = 0;
			for (int j = 0; j < fitFlux.size(); ++j) {
				rate += ionDamage.cascadeEfficiency[j] * fitFlux[j] *
						ionDamage.defectProduction[it][j];
			}
			rate = rate * fraction * fluxAmplitude / double(size);
			// Add it to the vector
			damageRate.push_back(rate);
		}
		damageRate.push_back(0.0);
		return damageRate;
	}

	/**
	 * Structure to hold damage dataFile
	 */
	struct ionDamageStruct {
		double perfectFraction;
		std::vector<double> fitRange;
	  std::vector<std::vector<double> > fitCoefficients;
		std::vector<double> cascadeEfficiency;
		std::vector<int> clusterSize;
		std::vector<std::vector<double> > defectProduction;
		std::vector<int> fluxIndex;
		std::vector<std::vector<double> > damageRate;
	} ionDamage;

public:

	/**
	 * The constructor
	 */
	AlloyFitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~AlloyFitFluxHandler() {}

	/**
	 * Compute and store the incident flux values at each grid point.
	 * \see IFluxHandler.h
	 */
	void initializeFluxHandler(IReactionNetwork *network,
			int surfacePos, std::vector<double> grid) {

		// Initialize the in-situ damage profile
		initializeInSituDamage();

		// Set the grid
		xGrid = grid;
		if (xGrid.size() == 0) return; // Treat this an an error for in-situ?

		// Iterate over all produced cluster species
		for (int it = 0; it < ionDamage.clusterSize.size(); ++it) {

			// Get size of cluster
			int size = ionDamage.clusterSize[it];

			// Check if cluster is interstitial type
			if (size > 0) {
				// See if theres an iType cluster of size
				auto fluxCluster = network->get(iType, size);
				if (fluxCluster) {
					fluxIndex = fluxCluster->getId() - 1;
					ionDamage.fluxIndex.push_back(fluxIndex);
					ionDamage.damageRate.push_back(AlloySetDamage(size,it,1.0));
				}
				// Otherwise the clusters must be frank and perfect type
				else {
					auto fluxCluster1 = network->get(frankType, size);
					auto fluxCluster2 = network->get(perfectType, size);
					if (!fluxCluster1 || !fluxCluster2) {
						// Throw error -> missing type
					}
					else {
						// Frank loop
						fluxIndex = fluxCluster1->getId() - 1;
						ionDamage.fluxIndex.push_back(fluxIndex);
						double frac = 1.0 - ionDamage.perfectFraction;
						ionDamage.damageRate.push_back(AlloySetDamage(size,it,frac));
						// Perfect loop
						fluxIndex = fluxCluster2->getId() - 1;
						ionDamage.fluxIndex.push_back(fluxIndex);
						frac = ionDamage.perfectFraction;
						ionDamage.damageRate.push_back(AlloySetDamage(size,it,frac));
					}
				}
			}
			// Check if cluster is vacancy type
			else if (size < 0) {
				size = -size;
				// See if theres an vType cluster of size
				auto fluxCluster = network->get(vType, size);
				if (fluxCluster) {
					fluxIndex = fluxCluster->getId() - 1;
					ionDamage.fluxIndex.push_back(fluxIndex);
					ionDamage.damageRate.push_back(AlloySetDamage(size,it,1.0));
				}
				// Otherwise the clusters must be faulted type
				else {
					fluxCluster = network->get(faultedType, size);
					if (!fluxCluster) {
						// Throw error -> no available type
					}
					else {
						// Faulted loop
						fluxIndex = fluxCluster->getId() - 1;
						ionDamage.fluxIndex.push_back(fluxIndex);
						ionDamage.damageRate.push_back(AlloySetDamage(size,it,1.0));
					}
				}
			}
			// Neither interstitial nor vacancy type cluster
			else {
				// Throw error for size 0 cluster
			}

		}

		// Check that the helium cluster is present in the network
		//if (!fluxCluster) {
		//	throw std::string(
		//			"\nThe single interstitial cluster is not present in the network, "
		//			"cannot use the flux option!");
		//}

		return;
	}

	/**
	 * This operation computes the flux due to incoming particles at a given grid point.
	 * \see IFluxHandler.h
	 */
	void computeIncidentFlux(double currentTime, double *updatedConcOffset, int xi, int surfacePos) {
		// Update the concentration array
		for (int it = 0; it < ionDamage.fluxIndex.size(); ++it) {
			updatedConcOffset[it] += ionDamage.damageRate[it][xi - surfacePos];
		}

		// There will be errors long before this check!
		//if (incidentFluxVec.size() == 0) {
		//	updatedConcOffset[fluxIndex] += fluxAmplitude;
		//	return;
		//}

		// Update the concentration array
		//updatedConcOffset[fluxIndex] += incidentFluxVec[xi - surfacePos];

		return;
	}

	void initializeInSituDamage() {
		ionDamage.perfectFraction = 0.8;
		ionDamage.fitRange = {0,	440};
		ionDamage.fitCoefficients =
		{
			{1.00835,	0.00524512,	-6.22E-05,	1.01E-06,	-6.17E-09,	1.44E-11,
					-1.14E-14},
			{1.15141,	0.00534056,	-5.59E-05,	1.16E-06,	-7.55E-09,	1.80E-11,
					-1.45E-14},
			{1.06211,	0.00581535,	-5.68E-05,	1.09E-06,	-7.06E-09,	1.67E-11,
					-1.34E-14},
			{9.70508,	0.00541518,	0.000425784,	-1.49E-06,	-9.90E-09,	4.47E-11,
					-4.60E-14}
		};
		ionDamage.cascadeEfficiency = {0.4,	0.25,	0.134,	0.12};
		ionDamage.clusterSize =
		{
			45,	40, 35,	30,	25,	20,	16,	12,	9,	8,	7,	6,	5,	4,	3,	2,	1,
			-1,	-2,	-3,	-4,	-5,	-9,	-12,	-15,	-20,	-25,	-30
		};
		ionDamage.defectProduction =
		{
			{0.000,		0.000,		0.013,		0.024},
			{0.000,		0.000,		0.015,		0.025},
			{0.000,		0.000,		0.016,		0.027},
			{0.000,		0.000,		0.018,		0.029},
			{0.000,		0.000,		0.021,		0.032},
			{0.000,		0.021,		0.024,		0.035},
			{0.000,		0.026,		0.029,		0.040},
			{0.000,		0.033,		0.036,		0.046},
			{0.000,		0.042,		0.044,		0.053},
			{0.040,		0.047,		0.049,		0.056},
			{0.054,		0.053,		0.054,		0.060},
			{0.063,		0.061,		0.060,		0.065},
			{0.076,		0.071,		0.069,		0.071},
    	{0.094,		0.086,		0.082,		0.079},
    	{0.124,		0.111,		0.101,		0.091},
    	{0.185,		0.158,		0.137,		0.112},
    	{0.364,		0.291,		0.232,		0.155},
			{0.673,		0.611,		0.537,		0.375},
    	{0.174,		0.181,		0.190,		0.187},
			{0.079,		0.089,		0.103,		0.125},
    	{0.045,		0.054,		0.067,		0.094},
    	{0.029,		0.036,		0.048,		0.075},
			{0.000,		0.013,		0.020,		0.042},
			{0.000,		0.008,		0.013,		0.031},
			{0.000,		0.005,		0.009,		0.025},
			{0.000,		0.003,		0.006,		0.019},
			{0.000,		0.000,		0.004,		0.015},
			{0.000,		0.000,		0.003,		0.012}
		};
		return;
	}

};
//end class AlloyFitFluxHandler

}

#endif
