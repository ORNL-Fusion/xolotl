#ifndef TEMPERATUREGRADIENTHANDLER_H
#define TEMPERATUREGRADIENTHANDLER_H

#include "ITemperatureHandler.h"

namespace xolotlCore{

/**
 * This class realizes the ITemperatureHandler, it is responsible for the
 * handling of a temperature constant with time but changing with location.
 */
class TemperatureGradientHandler: public ITemperatureHandler {

private:

	/**
	 * The surface temperature in Kelvin
	 */
	double surfaceTemperature;

	/**
	 * The temperature gradient
	 */
	double gradient;

	/**
	 * The default constructor is private because the TemperatureHandler
	 * must be initialized with a temperature
	 */
	TemperatureGradientHandler() :
		surfaceTemperature(0.0), gradient(0.0) {}

public:

	/**
	 * The constructor
	 *
	 * @param temp The surface temperature
	 * @param grad The temperature gradient
	 */
	TemperatureGradientHandler(double temp, double grad) :
		surfaceTemperature(temp), gradient(grad) {}

	/**
	 * The destructor.
	 */
	virtual ~TemperatureGradientHandler() {}

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
	virtual double getTemperature(const std::vector<double>& position,
			double) const {return surfaceTemperature - position[0] * gradient;}

}; //end class TemperatureGradientHandler

}

#endif
