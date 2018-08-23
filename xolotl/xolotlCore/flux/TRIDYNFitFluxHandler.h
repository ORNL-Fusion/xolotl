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
	 * Parameters for the Weibull fit that will be read from a file
	 */
	double A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15;

	/**
	 * The index of the V cluster.
	 */
	int vDefectIndex;

	/**
	 * The index of the I cluster.
	 */
	int iDefectIndex;

	/**
	 * Function that calculate the flux at a given position x (in nm).
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

public:

	/**
	 * The constructor
	 */
	TRIDYNFitFluxHandler() :
			A0(0.0), A1(0.0), A2(0.0), A3(0.0), A4(0.0), A5(0.0), A6(0.0), A7(
					0.0), A8(0.0), A9(0.0), A10(0.0), A11(0.0), A12(0.0), A13(
					0.0), A14(0.0), A15(0.0), vDefectIndex(-1), iDefectIndex(-1) {
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

			// Break the line into tokens
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
		}

		// Close the file
		paramFile.close();

		// Call the general method
		FluxHandler::initializeFluxHandler(network, surfacePos, grid);

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
//		// Set the I index corresponding the the single interstitial cluster here
//		auto iCluster = network->get(iType, 1);
//		// Check that the V cluster is present in the network
//		if (!iCluster) {
//			throw std::string(
//					"\nThe single interstitial cluster is not present in the network, "
//							"cannot use the flux option!");
//		}
//		iDefectIndex = iCluster->getId() - 1;

		return;
	}

	/**
	 * This operation computes the flux due to incoming particles at a given grid point.
	 * \see IFluxHandler.h
	 */
	void computeIncidentFlux(double currentTime, double *updatedConcOffset,
			int xi, int surfacePos) {
		// Recompute the flux vector if a time profile is used
		if (useTimeProfile) {
			fluxAmplitude = getProfileAmplitude(currentTime);
			recomputeFluxHandler(surfacePos);
		}

		// Update the concentration array
		updatedConcOffset[fluxIndex] += incidentFluxVec[xi - surfacePos];
//		updatedConcOffset[vDefectIndex] += incidentFluxVec[xi - surfacePos] * 4.25e-7;
//		updatedConcOffset[iDefectIndex] += incidentFluxVec[xi - surfacePos] * 4.25e-7;

		return;
	}

};
//end class TRIDYNFitFluxHandler

}

#endif
