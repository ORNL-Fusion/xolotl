#ifndef TRIDYNFITFLUXHANDLER_H
#define TRIDYNFITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>
#include <TokenizedLineReader.h>
#include <string.h>
#include <iostream>
#include <fstream>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for a (100) oriented tungsten material using TRIDYN input data.
 */
class TRIDYNFitFluxHandler: public FluxHandler {
private:

	/**
	 * Parameters for the polynomial fit that will be read from a file for He
	 */
	double A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15;

	/**
	 * Parameters for the polynomial fit that will be read from a file for W redeposition
	 */
	double B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15;

	/**
	 * The index of the V cluster.
	 */
	int vDefectIndex;

	/**
	 * The index of the I cluster.
	 */
	int iDefectIndex;

	/**
	 * The reduction factor between He and W redeposition.
	 */
	double reductionFactor;

	/**
	 * Value of the fit function integrated on the grid for W.
	 */
	double wNormFactor;

	/**
	 * Vector to hold the incident flux for W values at each grid
	 * point (x position).
	 */
	std::vector<double> incidentWFluxVec;

	/**
	 * Function that calculate the flux of He at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (100).
	 *
	 * @param x The position where to evaluate he fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		if (x > 11.5)
			return 0.0;

		// Compute the polynomial fit
		double value = A0 + A1 * x + A2 * pow(x, 2.0) + A3 * pow(x, 3.0)
				+ A4 * pow(x, 4.0) + A5 * pow(x, 5.0) + A6 * pow(x, 6.0)
				+ A7 * pow(x, 7.0) + A8 * pow(x, 8.0) + A9 * pow(x, 9.0)
				+ A10 * pow(x, 10.0) + A11 * pow(x, 11.0) + A12 * pow(x, 12.0)
				+ A13 * pow(x, 13.0) + A14 * pow(x, 14.0) + A15 * pow(x, 15.0);

		return std::max(value, 0.0);
	}

	/**
	 * Function that calculate the flux of W at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (100).
	 *
	 * @param x The position where to evaluate he fit
	 * @return The evaluated value
	 */
	double WFitFunction(double x) {
		// Compute the polynomial fit
		double value = B0 + B1 * x + B2 * pow(x, 2.0) + B3 * pow(x, 3.0)
				+ B4 * pow(x, 4.0) + B5 * pow(x, 5.0) + B6 * pow(x, 6.0)
				+ B7 * pow(x, 7.0) + B8 * pow(x, 8.0) + B9 * pow(x, 9.0)
				+ B10 * pow(x, 10.0) + B11 * pow(x, 11.0) + B12 * pow(x, 12.0)
				+ B13 * pow(x, 13.0) + B14 * pow(x, 14.0) + B15 * pow(x, 15.0);

		return std::max(value, 0.0);
	}

public:

	/**
	 * The constructor
	 */
	TRIDYNFitFluxHandler() : A0(0.0), A1(0.0), A2(0.0), A3(0.0), A4(0.0), A5(0.0), A6(0.0),
	A7(0.0), A8(0.0), A9(0.0), A10(0.0), A11(0.0), A12(0.0), A13(0.0), A14(0.0), A15(0.0), B0(0.0), B1(0.0), B2(0.0), B3(0.0), B4(0.0), B5(0.0), B6(0.0),
	B7(0.0), B8(0.0), B9(0.0), B10(0.0), B11(0.0), B12(0.0), B13(0.0), B14(0.0), B15(0.0),
	vDefectIndex(-1), iDefectIndex(-1), reductionFactor(0.0), wNormFactor(0.0) {
	}

	/**
	 * The Destructor
	 */
	~TRIDYNFitFluxHandler() {
	}

	/**
	 * Compute and store the incident flux values at each grid point.
	 * \see IFluxHandler.h
	 */
	void initializeFluxHandler(IReactionNetwork *network, int surfacePos,
			std::vector<double> grid) {
		// Read the parameter file
		std::ifstream paramFile;
		paramFile.open("tridyn.dat");

		if (!paramFile.good()) {
			std::cout
					<< "No parameter files for TRIDYN flux, the flux will be 0"
					<< std::endl;
		} else {
			// Get the line
			std::string line;
			getline(paramFile, line);

			// Build an input stream from the line string
			xolotlCore::TokenizedLineReader<double> reader;
			auto lineSS = std::make_shared<std::istringstream>(line);
			reader.setInputStream(lineSS);

			// Break the line into tokens for He
			auto tokens = reader.loadLine();
			// Set the parameters for the fit
			A0 = tokens[0];
			A1 = tokens[1];
			A2 = tokens[2];
			A3 = tokens[3];
			A4 = tokens[4];
			A5 = tokens[5];
			A6 = tokens[6];
			A7 = tokens[7];
			A8 = tokens[8];
			A9 = tokens[9];
			A10 = tokens[10];
			A11 = tokens[11];
			A12 = tokens[12];
			A13 = tokens[13];
			A14 = tokens[14];
			A15 = tokens[15];

			// Break the line into tokens for reduction factor
			getline(paramFile, line);
			lineSS = std::make_shared<std::istringstream>(line);
			reader.setInputStream(lineSS);
			tokens = reader.loadLine();
			reductionFactor = tokens[0];

			// Break the line into tokens for W redeposition
			getline(paramFile, line);
			lineSS = std::make_shared<std::istringstream>(line);
			reader.setInputStream(lineSS);
			tokens = reader.loadLine();
			// Set the parameters for the fit
			B0 = tokens[0];
			B1 = tokens[1];
			B2 = tokens[2];
			B3 = tokens[3];
			B4 = tokens[4];
			B5 = tokens[5];
			B6 = tokens[6];
			B7 = tokens[7];
			B8 = tokens[8];
			B9 = tokens[9];
			B10 = tokens[10];
			B11 = tokens[11];
			B12 = tokens[12];
			B13 = tokens[13];
			B14 = tokens[14];
			B15 = tokens[15];
		}

		// Close the file
		paramFile.close();

		// Call the general method
		FluxHandler::initializeFluxHandler(network, surfacePos, grid);


		// Compute the norm factor because the fit function has an
		// arbitrary amplitude
		wNormFactor = 0.0;
		// Loop on the x grid points skipping the first after the surface position
		// and last because of the boundary conditions
		for (int i = surfacePos + 1; i < xGrid.size() - 1; i++) {
			// Get the x position
			double x = xGrid[i] - xGrid[surfacePos];

			// Add the the value of the function times the step size
			wNormFactor += WFitFunction(x) * (xGrid[i] - xGrid[i-1]);
		}

		// Factor the incident flux will be multiplied by to get
		// the wanted intensity
		double fluxNormalized = 0.0;
		if (wNormFactor > 0.0) fluxNormalized = fluxAmplitude * reductionFactor / wNormFactor;

		// Clear the flux vector
		incidentWFluxVec.clear();
		// The first value corresponding to the surface position should always be 0.0
		incidentWFluxVec.push_back(0.0);

		// Starts a i = surfacePos + 1 because the first value was already put in the vector
		for (int i = surfacePos + 1; i < xGrid.size() - 1; i++) {
			// Get the x position
			auto x = xGrid[i] - xGrid[surfacePos];

			// Compute the flux value
			double incidentFlux = fluxNormalized * WFitFunction(x);
			// Add it to the vector
			incidentWFluxVec.push_back(incidentFlux);
		}

		// The last value should always be 0.0 because of boundary conditions
		incidentWFluxVec.push_back(0.0);

		// Set the flux index corresponding the the single helium cluster here
		auto fluxCluster = network->get(heType, 1);
		// Check that the helium cluster is present in the network
		if (!fluxCluster) {
			throw std::string(
					"\nThe single helium cluster is not present in the network, "
							"cannot use the flux option!");
		}
		fluxIndex = fluxCluster->getId() - 1;

//		// Set the V index corresponding the the single vacancy cluster here
//		auto vCluster = network->get(vType, 1);
//		// Check that the V cluster is present in the network
//		if (!vCluster) {
//			throw std::string(
//					"\nThe single vacancy cluster is not present in the network, "
//							"cannot use the flux option!");
//		}
//		vDefectIndex = vCluster->getId() - 1;
//
		// Set the I index corresponding the the single interstitial cluster here
		auto iCluster = network->get(iType, 1);
		// Check that the V cluster is present in the network
		if (!iCluster) {
			throw std::string(
					"\nThe single interstitial cluster is not present in the network, "
							"cannot use the flux option!");
		}
		iDefectIndex = iCluster->getId() - 1;

		return;
	}

	/**
	 * This operation computes the flux due to incoming particles at a given grid point.
	 * \see IFluxHandler.h
	 */
	void computeIncidentFlux(double currentTime, double *updatedConcOffset, int xi, int surfacePos) {
		// Recompute the flux vector if a time profile is used
		if (useTimeProfile) {
			fluxAmplitude = getProfileAmplitude(currentTime);
			recomputeFluxHandler(surfacePos);
		}

		// Update the concentration array
		updatedConcOffset[fluxIndex] += incidentFluxVec[xi - surfacePos];
//		updatedConcOffset[vDefectIndex] += incidentFluxVec[xi - surfacePos] * 4.25e-7;
		updatedConcOffset[iDefectIndex] += incidentWFluxVec[xi - surfacePos];

		return;
	}

};
//end class TRIDYNFitFluxHandler

}

#endif
