#ifndef A800H9MEVFLUXHANDLER_H
#define A800H9MEVFLUXHANDLER_H

#include <cmath>

#include <xolotl/core/flux/AlloyFluxHandler.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the AlloyFluxHandler interface for a 9 MeV ion flux.
 */
class A800H9MeVFluxHandler : public AlloyFluxHandler
{
public:
	/**
	 * The constructor
	 */
	A800H9MeVFluxHandler(const options::IOptions& options) :
		AlloyFluxHandler(options)
	{
		this->fluxI = {0, 0.960209, 0.00437383, 0.0017928, 0.00125702,
			0.00037196, 0.000410079, 0.000394621, 0.000148889, 0.000105098,
			0.000165166, 0, 0, 0.000113856, 0, 0, 0, 0, 0.000122614, 0, 0, 0, 0,
			0.000113856, 0, 0, 0, 0, 7.00654e-05, 0, 0, 0, 0, 1.31373e-05, 0, 0,
			0, 0, 1.31373e-05, 0, 0, 0, 0, 1.31373e-05};

		this->fluxV = {0, 0.980677, 0.00132315, 0.000405436, 0.000312638,
			0.000209156, 0.00010992, 0.000222715, 6.40205e-05, 6.76503e-05,
			6.77392e-05, 0, 0, 0.000112523, 0, 0, 0, 0, 2.72232e-05, 0, 0, 0, 0,
			5.44464e-05, 0, 0, 0, 0, 2.35934e-05, 0, 0, 0, 0, 3.26678e-05, 0, 0,
			0, 0, 1.63339e-05, 0, 0, 0, 0, 1.08893e-05, 0, 0, 0, 0, 2.54083e-05,
			0, 0, 0, 0, 1.4519e-05, 0, 0, 0, 0, 1.08893e-05, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 7.25952e-06};
	}

	/**
	 * The Destructor
	 */
	~A800H9MeVFluxHandler()
	{
	}
};
// end class A800H9MeVFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
