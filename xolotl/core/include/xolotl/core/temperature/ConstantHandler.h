#pragma once

#include <xolotl/core/temperature/TemperatureHandler.h>
#include <xolotl/options/IOptions.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
/**
 * This class realizes the TemperatureHandler, it is responsible for the
 * handling of a temperature constant with time and space.
 */
class ConstantHandler : public TemperatureHandler
{
private:
	/**
	 * The temperature in Kelvin
	 */
	double temperature;

public:
	/**
	 * Construct with provided constant temperature
	 *
	 * @param constTemperature the temperature
	 */
	ConstantHandler(double constTemperature);

	/**
	 * Construct with options
	 */
	ConstantHandler(const options::IOptions& options);

	/**
	 * The destructor.
	 */
	virtual ~ConstantHandler();

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 * Here it is a constant temperature.
	 *
	 * \see ITemperatureHandler.h
	 */
	double
	getTemperature(const plsm::SpaceVector<double, 3>&, double) const override
	{
		return temperature;
	}
};
} // namespace temperature
} // namespace core
} // namespace xolotl
