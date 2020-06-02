#ifndef TEMPERATUREPROFILEHANDLER_H
#define TEMPERATUREPROFILEHANDLER_H

#include <string>
#include <iostream>
#include <fstream>
#include <xolotl/core/temperature/ITemperatureHandler.h>

namespace xolotl {
namespace core {
namespace temperature {

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
	virtual void initializeTemperature(const int _dof,
			network::IReactionNetwork::SparseFillMap& ofillMap,
			network::IReactionNetwork::SparseFillMap& dfillMap) {

		// Set dof
		dof = _dof;

		// Add the temperature to ofill
		ofillMap[dof].emplace_back(dof);

		// Add the temperature to dfill
		dfillMap[dof].emplace_back(dof);

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
	virtual double getTemperature(const Point<3>&, double currentTime) const {
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
	 * This operation sets the heat conductivity to use in the equation.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void setHeatConductivity(double cond) {
		return;
	}

	/**
	 * This operation sets the surface position.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void updateSurfacePosition(int surfacePos) {
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
			double *updatedConcOffset, double hxLeft, double hxRight, int xi,
			double sy = 0.0, int iy = 0, double sz = 0.0, int iz = 0) {
		return;
	}

	/**
	 * Compute the partials due to the heat equation.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual bool computePartialsForTemperature(double *val, int *indices,
			double hxLeft, double hxRight, int xi, double sy = 0.0, int iy = 0,
			double sz = 0.0, int iz = 0) {
		return false;
	}

};
//end class TemperatureProfileHandler

}
}
}

#endif
