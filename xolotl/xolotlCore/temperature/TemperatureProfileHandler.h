#ifndef TEMPERATUREPROFILEHANDLER_H
#define TEMPERATUREPROFILEHANDLER_H

#include "ITemperatureHandler.h"
#include <string>
#include <iostream>
#include <fstream>

namespace xolotlCore{

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
	 * The default constructor is private because the TemperatureProfileHandler
	 * must be initialized with an input file.
	 */
	TemperatureProfileHandler() :
		tempFile("") {}

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
	TemperatureProfileHandler(std::string profileFileName)
		: tempFile(profileFileName) {}

	/**
	 * The destructor.
	 */
	virtual ~TemperatureProfileHandler() {}

	/**
	 * This operation reads in the time and temperature data from the input
	 * temperature file that was specified by the command line
	 */
	virtual void initializeTemperature() {
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
	 * This operation linearly interpolates the data read from the input
	 * temperature file and returns the temperature at the given position
	 * and time.
	 *
	 * @param position The position
	 * @param currentTime The time
	 * @return The temperature
	 */
	virtual double getTemperature(std::vector<double> position, double currentTime) const {
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
		for (int k = 0; k < time.size() - 1; k++) {
			if (currentTime < time[k]) continue;
			if (currentTime > time[k + 1]) continue;

			// Compute the amplitude following a linear interpolation between
			// the two stored values
			f = temp[k]
					+ (temp[k + 1] - temp[k]) * (currentTime - time[k])
							/ (time[k + 1] - time[k]);
			break;
		}

		return f;
	}

}; //end class TemperatureProfileHandler

}

#endif
