#ifndef TEMPERATUREHANDLER_H
#define TEMPERATUREHANDLER_H

#include "ITemperatureHandler.h"

namespace xolotlSolver{

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * and outgoing Temperature calculations.
 */
class TemperatureHandler: public ITemperatureHandler {

public:

	/**
	 * The constructor
	 */
	TemperatureHandler();

	/**
	 * The Destructor
	 */
	~TemperatureHandler();

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 * @param position        The position
	 * @param currentTime     The time
	 * @return temperature   The temperature
	 */
	virtual double getTemperature(std::vector<double> position, double currentTime) const;

}; //end class TemperatureHandler

}

#endif
