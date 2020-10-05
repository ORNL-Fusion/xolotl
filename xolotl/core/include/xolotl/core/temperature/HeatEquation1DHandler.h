#ifndef HEATEQUATION1DHANDLER_H
#define HEATEQUATION1DHANDLER_H

#include <xolotl/core/temperature/HeatEquationHandler.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
/**
 * This class realizes the ITemperatureHandler, it is responsible for the
 * handling of the heat equation in 1D.
 */
class HeatEquation1DHandler : public HeatEquationHandler
{
public:
	HeatEquation1DHandler() = delete;

	/**
	 * The constructor
	 *
	 * @param flux The heat flux
	 * @param bulkTemp The temperature in the bulk
	 */
	HeatEquation1DHandler(double flux, double bulkTemp) :
		HeatEquationHandler(flux, bulkTemp, 1)
	{
	}

	/**
	 * The destructor.
	 */
	virtual ~HeatEquation1DHandler()
	{
	}
};
// end class HeatEquation1DHandler

} // namespace temperature
} // namespace core
} // namespace xolotl

#endif
