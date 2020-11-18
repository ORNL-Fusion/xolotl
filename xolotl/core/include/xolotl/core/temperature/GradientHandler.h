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
 * handling of a temperature constant with time but changing with location.
 */
class GradientHandler : public TemperatureHandler
{
private:
	/**
	 * The surface temperature in Kelvin
	 */
	double surfaceTemperature;

	/**
	 * The bulk temperature in Kelvin
	 */
	double bulkTemperature;

public:
	GradientHandler() = delete;

	/**
	 * The constructor
	 *
	 * @param temp The surface temperature
	 * @param grad The bulk temperature
	 */
	GradientHandler(double surfTemp, double bulkTemp);

	/**
	 * Construct from options
	 */
	GradientHandler(const options::IOptions& options);

	/**
	 * The destructor.
	 */
	virtual ~GradientHandler();

	/**
	 * \see ITemperatureHandler.h
	 */
	double
	getTemperature(
		const plsm::SpaceVector<double, 3>& fraction, double) const override;
};
// end class GradientHandler

} // namespace temperature
} // namespace core
} // namespace xolotl
