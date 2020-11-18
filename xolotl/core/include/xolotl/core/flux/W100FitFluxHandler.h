#ifndef W100FITFLUXHANDLER_H
#define W100FITFLUXHANDLER_H

#include <cmath>

#include <xolotl/core/flux/PSIFluxHandler.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the PSIFluxHandler interface to calculate the incident
 * helium flux for a (100) oriented tungsten material.
 */
class W100FitFluxHandler : public PSIFluxHandler
{
private:
	/**
	 * \see FluxHandler.h
	 */
	double
	FitFunction(double x)
	{
		// Value at which the flux goes to 0
		double x1 = 10.0;

		if (x > x1)
			return 0.0;

		// Compute the fit
		double value = 7.00876507 + 0.6052078 * x - 3.01711048 * pow(x, 2) +
			1.36595786 * pow(x, 3) - 0.295595 * pow(x, 4) +
			0.03597462 * pow(x, 5) - 0.0025142 * pow(x, 6) +
			0.0000942235 * pow(x, 7) - 0.0000014679 * pow(x, 8);

		return value;
	}

public:
	/**
	 * The constructor
	 */
	W100FitFluxHandler(const options::IOptions& options) :
		PSIFluxHandler(options)
	{
	}

	/**
	 * The Destructor
	 */
	~W100FitFluxHandler()
	{
	}
};
// end class W100FitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
