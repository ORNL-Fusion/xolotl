#ifndef TRIDYNFITFLUXHANDLER_H
#define TRIDYNFITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>
#include <TokenizedLineReader.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <MPIUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for a (100) oriented tungsten material using TRIDYN input data.
 */
class TRIDYNFitFluxHandler: public FluxHandler {
private:

	/**
	 * Parameters for the polynomial fits that will be read from a file
	 */
	std::vector<std::vector<double> > fitParams;

	/**
	 * The total depths on which the flux is defined.
	 */
	std::vector<double> totalDepths;

	/**
	 * The reduction factor between He and other depositions.
	 */
	std::vector<double> reductionFactors;

	/**
	 * Value of the fit function integrated on the grid.
	 */
	std::vector<double> normFactors;

	/**
	 * Vector to hold the incident flux for W values at each grid
	 * point (x position).
	 */
	std::vector<double> incidentWFluxVec;

	/**
	 * Vector to hold the incident flux for D values at each grid
	 * point (x position).
	 */
	std::vector<double> incidentDFluxVec;

	/**
	 * Vector to hold the incident flux for T values at each grid
	 * point (x position).
	 */
	std::vector<double> incidentTFluxVec;

	/**
	 * Function that calculate the flux of He at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (100).
	 *
	 * @param x The position where to evaluate the fit
	 * @param i The indice of the fit
	 * @return The evaluated value
	 */
	double FitFunction(double x, int i) {
		if (x > totalDepths[i])
			return 0.0;

		// Compute the polynomial fit
		auto params = fitParams[i];
		double value = params[0] + params[1] * x + params[2] * pow(x, 2.0)
				+ params[3] * pow(x, 3.0) + params[4] * pow(x, 4.0)
				+ params[5] * pow(x, 5.0) + params[6] * pow(x, 6.0)
				+ params[7] * pow(x, 7.0) + params[8] * pow(x, 8.0)
				+ params[9] * pow(x, 9.0) + params[10] * pow(x, 10.0)
				+ params[11] * pow(x, 11.0) + params[12] * pow(x, 12.0)
				+ params[13] * pow(x, 13.0) + params[14] * pow(x, 14.0)
				+ params[15] * pow(x, 15.0);

		return std::max(value, 0.0);
	}

public:

	/**
	 * The constructor
	 */
	TRIDYNFitFluxHandler() {
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
	void initializeFluxHandler(const IReactionNetwork& network, int surfacePos,
			std::vector<double> grid) {
		// Read the parameter file
		std::ifstream paramFile;
		paramFile.open("tridyn.dat");

		if (!paramFile.good()) {
			// Print a message
			std::cout
					<< "No parameter files for TRIDYN flux, the flux will be 0"
					<< std::endl;

			// Set the depths to 0.0
			totalDepths.push_back(0.0);
			totalDepths.push_back(0.0);
			totalDepths.push_back(0.0);
			totalDepths.push_back(0.0);

			// Set the reduction factors to 0.0
			reductionFactors.push_back(0.0);
			reductionFactors.push_back(0.0);
			reductionFactors.push_back(0.0);
			reductionFactors.push_back(0.0);
		} else {
			// Build an input stream from the string
			xolotlCore::TokenizedLineReader<double> reader;
			// Get the line
			std::string line;
			getline(paramFile, line);
			auto lineSS = std::make_shared<std::istringstream>(line);
			reader.setInputStream(lineSS);

			// Read the first line
			auto tokens = reader.loadLine();
			// And start looping on the lines
			int i = 0;
			while (i < 4) {
				// Get the fraction
				reductionFactors.push_back(tokens[0]);

				// Set the parameters for the fit
				getline(paramFile, line);
				lineSS = std::make_shared<std::istringstream>(line);
				reader.setInputStream(lineSS);
				tokens = reader.loadLine();
				std::vector<double> params;
				params.push_back(tokens[0]);
				params.push_back(tokens[1]);
				params.push_back(tokens[2]);
				params.push_back(tokens[3]);
				params.push_back(tokens[4]);
				params.push_back(tokens[5]);
				params.push_back(tokens[6]);
				params.push_back(tokens[7]);
				params.push_back(tokens[8]);
				params.push_back(tokens[9]);
				params.push_back(tokens[10]);
				params.push_back(tokens[11]);
				params.push_back(tokens[12]);
				params.push_back(tokens[13]);
				params.push_back(tokens[14]);
				params.push_back(tokens[15]);
				fitParams.push_back(params);
				// Set the total depth where the fit is defined
				totalDepths.push_back(tokens[16] + 0.1);

				// Read the next line
				getline(paramFile, line);
				lineSS = std::make_shared<std::istringstream>(line);
				reader.setInputStream(lineSS);
				tokens = reader.loadLine();
				// Increase the loop number
				i++;
			}
		}

		// Close the file
		paramFile.close();

		// Set the grid
		xGrid = grid;

		if (xGrid.size() == 0)
			return;

		// Clear the norm factors
		normFactors.clear();

		// Loop on the different types of clusters, He, W, D, T
		for (int index = 0; index < 4; index++) {
			// Compute the norm factor because the fit function has an
			// arbitrary amplitude
			normFactors.push_back(0.0);
			// Loop on the x grid points skipping the first after the surface position
			// and last because of the boundary conditions
			for (int i = surfacePos + 1; i < xGrid.size() - 3; i++) {
				// Get the x position
				double x = xGrid[i + 1] - xGrid[surfacePos + 1];

				// Add the the value of the function times the step size
				normFactors[index] += FitFunction(x, index)
						* (xGrid[i + 1] - xGrid[i]);
			}

			// Factor the incident flux will be multiplied by to get
			// the wanted intensity
			double fluxNormalized = 0.0;
			if (normFactors[index] > 0.0) {
				fluxNormalized = fluxAmplitude * reductionFactors[index]
						/ normFactors[index];
			}

			// Select the right vector
			std::vector<double> *fluxVec = nullptr;
			switch (index) {
			case 0:
				fluxVec = &incidentFluxVec;
				break;
			case 1:
				fluxVec = &incidentWFluxVec;
				break;
			case 2:
				fluxVec = &incidentDFluxVec;
				break;
			case 3:
				fluxVec = &incidentTFluxVec;
				break;
			default:
				break;
			}

			// Clear the flux vector
			fluxVec->clear();
			// The first value corresponding to the surface position should always be 0.0
			fluxVec->push_back(0.0);

			// Starts at i = surfacePos + 1 because the first value was already put in the vector
			for (int i = surfacePos + 1; i < xGrid.size() - 3; i++) {
				// Get the x position
				auto x = xGrid[i + 1] - xGrid[surfacePos + 1];

				// Compute the flux value
				double incidentFlux = fluxNormalized * FitFunction(x, index);
				// Add it to the vector
				fluxVec->push_back(incidentFlux);
			}

			// The last value should always be 0.0 because of boundary conditions
			fluxVec->push_back(0.0);
		}

		// Set the flux index corresponding the the single helium cluster here
		auto fluxCluster = network.get(Species::He, 1);
		// Check that the helium cluster is present in the network
		if (!fluxCluster) {
			// Set its flux to 0.0
			fluxIndices.push_back(0);
			if (reductionFactors[0] > 0.0) {
				throw std::string(
						"\nThe single helium cluster is not present in the network, "
								"cannot use the flux option!");
			}
		} else
			fluxIndices.push_back(fluxCluster->getId() - 1);
//
		// Set the I index corresponding the the single interstitial cluster here
		fluxCluster = network.get(Species::I, 1);
		// Check that the I cluster is present in the network
		if (!fluxCluster) {
			// Set its flux to 0.0
			fluxIndices.push_back(0);
			if (reductionFactors[1] > 0.0) {
				throw std::string(
						"\nThe single interstitial cluster is not present in the network, "
								"cannot use the flux option!");
			}
		} else
			fluxIndices.push_back(fluxCluster->getId() - 1);

		// Set the D index corresponding the the single deuterium cluster here
		fluxCluster = network.get(Species::D, 1);
		// Check that the D cluster is present in the network
		if (!fluxCluster) {
			// Set its flux to 0.0
			fluxIndices.push_back(0);
			if (reductionFactors[2] > 0.0) {
				throw std::string(
						"\nThe single deuterium cluster is not present in the network, "
								"cannot use the flux option!");
			}
		} else
			fluxIndices.push_back(fluxCluster->getId() - 1);

		// Set the T index corresponding the the single tritium cluster here
		fluxCluster = network.get(Species::T, 1);
		// Check that the T cluster is present in the network
		if (!fluxCluster) {
			// Set its flux to 0.0
			fluxIndices.push_back(0);
			if (reductionFactors[3] > 0.0) {
				throw std::string(
						"\nThe single tritium cluster is not present in the network, "
								"cannot use the flux option!");
			}
		} else
			fluxIndices.push_back(fluxCluster->getId() - 1);

		// Prints both incident vectors in a file
		int procId;
		auto xolotlComm = MPIUtils::getMPIComm();
		MPI_Comm_rank(xolotlComm, &procId);
		if (procId == 0) {
			std::ofstream outputFile;
			outputFile.open("incidentVectors.txt");
			for (int i = 0; i < incidentFluxVec.size(); i++) {
				outputFile << grid[surfacePos + i + 1] - grid[surfacePos + 1]
						<< " " << incidentFluxVec[i] << " "
						<< incidentWFluxVec[i] << " " << incidentDFluxVec[i]
						<< " " << incidentTFluxVec[i] << std::endl;
			}
			outputFile.close();
		}

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
		updatedConcOffset[fluxIndices[0]] += incidentFluxVec[xi - surfacePos]; // He
//		updatedConcOffset[fluxIndices[1]] += incidentFluxVec[xi - surfacePos] * 4.25e-7;
		updatedConcOffset[fluxIndices[1]] += incidentWFluxVec[xi - surfacePos]; // I
		updatedConcOffset[fluxIndices[2]] += incidentDFluxVec[xi - surfacePos]; // D
		updatedConcOffset[fluxIndices[3]] += incidentTFluxVec[xi - surfacePos]; // T

		return;
	}

	void recomputeFluxHandler(int surfacePos) {
		// Loop on the different types of clusters, He, W, D, t
		for (int index = 0; index < 4; index++) {
			// Select the right vector
			std::vector<double> *fluxVec = nullptr;
			switch (index) {
			case 0:
				fluxVec = &incidentFluxVec;
				break;
			case 1:
				fluxVec = &incidentWFluxVec;
				break;
			case 2:
				fluxVec = &incidentDFluxVec;
				break;
			case 3:
				fluxVec = &incidentTFluxVec;
				break;
			default:
				break;
			}

			// Factor the incident flux will be multiplied by
			double fluxNormalized = 0.0;
			if (normFactors[index] > 0.0)
				fluxNormalized = fluxAmplitude * reductionFactors[index]
						/ normFactors[index];

			// Starts a i = surfacePos + 1 because the first values were already put in the vector
			for (int i = surfacePos + 1; i < xGrid.size() - 3; i++) {
				// Get the x position
				auto x = xGrid[i + 1] - xGrid[surfacePos + 1];

				// Compute the flux value
				double incidentFlux = fluxNormalized * FitFunction(x, index);
				// Add it to the vector
				fluxVec->at(i - surfacePos) = incidentFlux;
			}
		}

		return;
	}

};
//end class TRIDYNFitFluxHandler

}

#endif
