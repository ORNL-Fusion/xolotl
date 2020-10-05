#ifndef HEATEQUATION2DHANDLER_H
#define HEATEQUATION2DHANDLER_H

#include <xolotl/core/temperature/HeatEquationHandler.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
/**
 * This class realizes the ITemperatureHandler, it is responsible for the
 * handling of the heat equation in 2D.
 */
class HeatEquation2DHandler : public HeatEquationHandler
{
public:
	HeatEquation2DHandler() = delete;

	/**
	 * The constructor
	 *
	 * @param flux The heat flux
	 * @param bulkTemp The temperature in the bulk
	 */
	HeatEquation2DHandler(double flux, double bulkTemp) :
		HeatEquationHandler(flux, bulkTemp, 2)
	{
	}

	/**
	 * The destructor.
	 */
	virtual ~HeatEquation2DHandler()
	{
	}
};
// end class HeatEquation2DHandler

} // namespace temperature
} // namespace core
} // namespace xolotl

#endif
