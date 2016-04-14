#ifndef TEMPERATUREHANDLER_H
#define TEMPERATUREHANDLER_H

#include "ITemperatureHandler.h"

namespace xolotlCore{

/**
 * This class realizes the ITemperatureHandler, it is responsible for the
 * handling of a temperature constant with time.
 */
class TemperatureHandler: public ITemperatureHandler {

private:

	/**
	 * The temperature in Kelvin
	 */
	double temperature;

	/**
	 * The default constructor is private because the TemperatureHandler
	 * must be initialized with a temperature
	 */
	TemperatureHandler() :
		temperature(0.0) {}

public:

	/**
	 * The constructor
	 *
	 * @param constTemperature the temperature
	 */
	TemperatureHandler(double constTemperature) :
		temperature(constTemperature) {}

	/**
	 * The destructor.
	 */
	virtual ~TemperatureHandler() {}

	/**
	 * This operation reads in the time and temperature data from the input
	 * temperature file that was specified by the command line if the profile
	 * option is used;
	 */
	virtual void initializeTemperature() {}

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 *
	 * @return The temperature
	 */
	virtual double getTemperature(const std::vector<double>&,
			double) const {return temperature;}

}; //end class TemperatureHandler

}

#endif
