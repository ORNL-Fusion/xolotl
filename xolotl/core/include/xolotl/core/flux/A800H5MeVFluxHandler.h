#ifndef A800H5MEVFLUXHANDLER_H
#define A800H5MEVFLUXHANDLER_H

#include <cmath>

#include <xolotl/core/flux/AlloyFluxHandler.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the AlloyFluxHandler interface for a 5 MeV ion flux.
 */
class A800H5MeVFluxHandler : public AlloyFluxHandler
{
public:
	/**
	 * The constructor
	 */
	A800H5MeVFluxHandler(const options::IOptions& options) :
		AlloyFluxHandler(options)
	{
		this->fluxI = {0, 0.950727, 0.00534434, 0.00216457, 0.00151063,
			0.000464616, 0.000502879, 0.000475748, 0.00019176, 0.00013536,
			0.000204779, 0, 0, 0.00014664, 0, 0, 0, 0, 0.00015792, 0, 0, 0, 0,
			0.00014664, 0, 0, 0, 0, 9.02398e-05, 0, 0, 0, 0, 1.692e-05, 0, 0, 0,
			0, 1.692e-05, 0, 0, 0, 0, 1.692e-05};

		this->fluxV = {0, 0.975816, 0.00162583, 0.000492815, 0.000383084,
			0.000258614, 0.000136676, 0.000263355, 7.83439e-05, 8.30188e-05,
			8.09802e-05, 0, 0, 0.000144922, 0, 0, 0, 0, 3.50617e-05, 0, 0, 0, 0,
			7.01235e-05, 0, 0, 0, 0, 3.03868e-05, 0, 0, 0, 0, 4.20741e-05, 0, 0,
			0, 0, 2.1037e-05, 0, 0, 0, 0, 1.40247e-05, 0, 0, 0, 0, 3.27243e-05,
			0, 0, 0, 0, 1.86996e-05, 0, 0, 0, 0, 1.40247e-05, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 9.3498e-06};
	}

	/**
	 * The Destructor
	 */
	~A800H5MeVFluxHandler()
	{
	}
};
// end class A800H5MeVFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
