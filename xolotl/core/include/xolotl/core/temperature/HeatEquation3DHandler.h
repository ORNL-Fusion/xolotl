#ifndef HEATEQUATION3DHANDLER_H
#define HEATEQUATION3DHANDLER_H

#include <xolotl/core/temperature/HeatEquationHandler.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
/**
 * This class realizes the ITemperatureHandler, it is responsible for the
 * handling of the heat equation in 3D.
 */
class HeatEquation3DHandler : public HeatEquationHandler
{
public:
	HeatEquation3DHandler() = delete;

	/**
	 * The constructor
	 *
	 * @param flux The heat flux
	 * @param bulkTemp The temperature in the bulk
	 */
	HeatEquation3DHandler(double flux, double bulkTemp) :
		HeatEquationHandler(flux, bulkTemp, 3)
	{
	}

	/**
	 * The destructor.
	 */
	virtual ~HeatEquation3DHandler()
	{
	}
};
// end class HeatEquation3DHandler

} // namespace temperature
} // namespace core
} // namespace xolotl

#endif
