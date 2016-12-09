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
	double A;
	double L;
	double K;

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (100).
	 *
	 * @param x The position where to evaluate he fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		// Compute the Weibull fit
		double value = A * (K / L)
				* (pow((x / L), (K - 1.0)))
				* exp(-pow((x / L), K));

		return value;
	}

public:

	/**
	 * The constructor
	 */
	TRIDYNFitFluxHandler() : A(0.0), L(0.0), K(0.0) {
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
		A = tokens[0];
		L = tokens[1];
		K = tokens[2];

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

		return;
	}

};
//end class TRIDYNFitFluxHandler

}

#endif
