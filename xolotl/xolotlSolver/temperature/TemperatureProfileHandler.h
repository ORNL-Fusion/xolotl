#ifndef TEMPERATUREPROFILEHANDLER_H
#define TEMPERATUREPROFILEHANDLER_H

#include "ITemperatureHandler.h"
#include <string>

namespace xolotlSolver{

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * and outgoing Temperature calculations.
 */
class TemperatureProfileHandler: public ITemperatureHandler {

private:

	std::string tempFile;

	/**
	 * The default constructor is private because the TemperatureProfileHandler
	 * must be initialized with an input file
	 */
	TemperatureProfileHandler() :
		tempFile("")
	{ }

protected:

	/**
	 * Vector to hold the time read from the input
	 * temperature file
	 */
	std::vector<double> time;

	/**
	 * Vector to hold the temperature read from the input
	 * temperature file
	 */
	std::vector<double> temp;

public:

	/**
	 * The constructor
	 */
	TemperatureProfileHandler(std::string profileFileName)
		: tempFile(profileFileName)
	{ }

	/**
	 * The Destructor
	 */
	virtual ~TemperatureProfileHandler() { }

	/**
	 * This operation reads in the time and temperature data from the input
	 * temperature file that was specified by the command line
	 */
	virtual void initializeTemperature();

	/**
	 * This operation linearly interpolates the data read from the input
	 * temperature file and returns the temperature at the given position
	 * and time.
	 * @param position        The position
	 * @param currentTime     The time
	 * @return temperature   The temperature
	 */
	virtual double getTemperature(std::vector<double> position, double currentTime) const;

}; //end class TemperatureProfileHandler

}

#endif
