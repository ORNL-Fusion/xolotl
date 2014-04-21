#ifndef ITEMPERATUREHANDLER_H
#define ITEMPERATUREHANDLER_H

#include <vector>

namespace xolotlSolver{

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * and outgoing Temperature calculations.
 */
class ITemperatureHandler {

public:

	virtual ~ITemperatureHandler() { }

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 * @param position        The position
	 * @param currentTime     The time
	 * @return temperature   The temperature
	 */
	virtual double getTemperature(std::vector<double> position, double currentTime) const = 0;

}; //end class ITemperatureHandler

}

#endif
