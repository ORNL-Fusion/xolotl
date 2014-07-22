#ifndef TEMPERATUREHANDLER_H
#define TEMPERATUREHANDLER_H

#include "ITemperatureHandler.h"

namespace xolotlSolver{

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * and outgoing Temperature calculations.
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
		temperature( 0.0e-16 )
	{ }

public:

	/**
	 * The constructor
	 * @param constTemperature the temperature
	 */
	TemperatureHandler(double constTemperature) :
		temperature( constTemperature )
	{ }

	/**
	 * The Destructor
	 */
	virtual ~TemperatureHandler() { }

	virtual void initializeTemperature() { }

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 * @param position        The position
	 * @param currentTime     The time
	 * @return temperature   The temperature
	 */
	virtual double getTemperature(std::vector<double> position,
			double currentTime) const	{ return temperature; }

}; //end class TemperatureHandler

}

#endif
