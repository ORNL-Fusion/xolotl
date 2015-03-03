#ifndef ITEMPERATUREHANDLER_H
#define ITEMPERATUREHANDLER_H

#include <vector>
#include <memory>

namespace xolotlCore{

/**
 * Realizations of this interface are responsible for handling the temperature.
 */
class ITemperatureHandler {

public:

	/**
	 * The destructor.
	 */
	virtual ~ITemperatureHandler() { }

	/**
	 * This operation reads in the time and temperature data from the input
	 * temperature file that was specified by the command line if the profile
	 * option is used;
	 */
	virtual void initializeTemperature() = 0;

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 *
	 * @param position The position
	 * @param currentTime The time
	 * @return The temperature
	 */
	virtual double getTemperature(std::vector<double> position, double currentTime) const = 0;

}; //end class ITemperatureHandler

}

#endif
