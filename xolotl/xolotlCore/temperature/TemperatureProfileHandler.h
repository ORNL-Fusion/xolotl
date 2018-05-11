#ifndef TEMPERATUREPROFILEHANDLER_H
#define TEMPERATUREPROFILEHANDLER_H

#include "ITemperatureHandler.h"
#include <string>
#include <iostream>
#include <fstream>

namespace xolotlCore {

/**
 * This class realizes the ITemperatureHandler, it is responsible for the
 * handling of a temperature changing with time.
 */
class TemperatureProfileHandler: public ITemperatureHandler {

private:

	/**
	 * The name of the file were the profile is stored.
	 */
	std::string tempFile;

	/**
	 * The number of degrees of freedom in the network
	 */
	int dof;

	/**
	 * The default constructor is private because the TemperatureProfileHandler
	 * must be initialized with an input file.
	 */
	TemperatureProfileHandler() :
			tempFile(""), dof(0) {
	}

	/**
	 * Vector to hold the time read from the input
	 * temperature file.
	 */
	std::vector<double> time;

	/**
	 * Vector to hold the temperature read from the input
	 * temperature file.
	 */
	std::vector<double> temp;

public:

	/**
	 * The constructor.
	 *
	 * @param profileFileName The name of the profile file
	 */
	TemperatureProfileHandler(const std::string& profileFileName) :
			tempFile(profileFileName), dof(0) {
	}

	/**
	 * The destructor.
	 */
	virtual ~TemperatureProfileHandler() {
	}

	/**
	 * This operation initializes the ofill and dfill arrays so that the
	 * temperature is connected correctly in the solver.
	 * It also reads in the time and temperature data from the input
	 * temperature file that was specified by the command line.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void initializeTemperature(const IReactionNetwork& network,
            IReactionNetwork::SparseFillMap& ofillMap,
            IReactionNetwork::SparseFillMap& dfillMap) {

		// Set dof
		dof = network.getDOF();

		// Add the temperature to ofill
        ofillMap[(dof - 1)].emplace_back(dof - 1);

		// Add the temperature to dfill
        dfillMap[(dof - 1)].emplace_back(dof - 1);

		// Open file dataFile.dat containing the time and temperature
		std::ifstream inputFile(tempFile.c_str());
		std::string line;

		// Read the file and store the values in the two vectors
		while (getline(inputFile, line)) {
			if (!line.length() || line[0] == '#')
				continue;
			double xtemp = 0.0, ytemp = 0.0;
			sscanf(line.c_str(), "%lf %lf", &xtemp, &ytemp);
			time.push_back(xtemp);
			temp.push_back(ytemp);
		}

		return;
	}

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 * It linearly interpolates the data read from the input
	 * temperature file.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual double getTemperature(const Point<3>& position,
			double currentTime) const {
		// Initialize the value to return
		double f = 0.0;

		// If the time is smaller than or equal than the first stored time
		if (currentTime <= time[0])
			return f = temp[0];

		// If the time is larger or equal to the last stored time
		if (currentTime >= time[time.size() - 1])
			return f = temp[time.size() - 1];

		// Else loop to determine the interval the time falls in
		// i.e. time[k] < time < time[k + 1]
		for (unsigned int k = 0; k < time.size() - 1; k++) {
			if (currentTime < time[k])
				continue;
			if (currentTime > time[k + 1])
				continue;

			// Compute the amplitude following a linear interpolation between
			// the two stored values
			f = temp[k]
					+ (temp[k + 1] - temp[k]) * (currentTime - time[k])
							/ (time[k + 1] - time[k]);
			break;
		}

		return f;
	}

	/**
	 * This operation sets the temperature given by the solver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void setTemperature(double * solution) {
		return;
	}

	/**
	 * This operation sets the heat coefficient to use in the equation.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void setHeatCoefficient(double coef) {
		return;
	}

	/**
	 * This operation sets the surface position.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void updateSurfacePosition(double surfacePos) {
		return;
	}

	/**
	 * Compute the flux due to the heat equation.
	 * This method is called by the RHSFunction from the PetscSolver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void computeTemperature(double **concVector,
			double *updatedConcOffset, double hxLeft, double hxRight) {
		return;
	}

	/**
	 * Compute the partials due to the heat equation.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void computePartialsForTemperature(double *val, int *indices,
			double hxLeft, double hxRight) {
		// Set the cluster index, the PetscSolver will use it to compute
		// the row and column indices for the Jacobian
		indices[0] = dof - 1;

		// Compute the partial derivatives for diffusion of this cluster
		// for the middle, left, and right grid point
		val[0] = 0.0; // middle
		val[1] = 0.0; // left
		val[2] = 0.0; // right

		return;
	}

};
//end class TemperatureProfileHandler

}

#endif
