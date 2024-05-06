#ifndef A800HNEUTRONFLUXHANDLER_H
#define A800HNEUTRONFLUXHANDLER_H

#include <cmath>

#include <xolotl/core/flux/AlloyFluxHandler.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the AlloyFluxHandler interface for a neutron flux.
 */
class A800HNeutronFluxHandler : public AlloyFluxHandler
{
public:
	/**
	 * The constructor
	 */
	A800HNeutronFluxHandler(const options::IOptions& options) :
		AlloyFluxHandler(options)
	{
		this->fluxI = {0.0, 4.14E-01, 5.80E-02, 2.14E-02, 1.44E-02, 5.84E-03,
			5.60E-03, 4.65E-03, 2.86E-03, 2.02E-03, 2.46E-03, 0.0, 0.0,
			2.19E-03, 0.0, 0.0, 0.0, 0.0, 2.36E-03, 0.0, 0.0, 0.0, 0.0,
			2.19E-03, 0.0, 0.0, 0.0, 0.0, 1.35E-03, 0.0, 0.0, 0.0, 0.0,
			2.52E-04, 0.0, 0.0, 0.0, 0.0, 2.52E-04, 0.0, 0.0, 0.0, 0.0,
			2.52E-04};

		this->fluxV = {0.0, 6.93E-01, 1.83E-02, 5.13E-03, 4.23E-03, 3.04E-03,
			1.67E-03, 2.15E-03, 8.57E-04, 9.27E-04, 7.33E-04, 0.0, 0.0,
			2.16E-03, 0.0, 0.0, 0.0, 0.0, 5.23E-04, 0.0, 0.0, 0.0, 0.0,
			1.05E-03, 0.0, 0.0, 0.0, 0.0, 4.53E-04, 0.0, 0.0, 0.0, 0.0,
			6.28E-04, 0.0, 0.0, 0.0, 0.0, 3.14E-04, 0.0, 0.0, 0.0, 0.0,
			2.09E-04, 0.0, 0.0, 0.0, 0.0, 4.88E-04, 0.0, 0.0, 0.0, 0.0,
			2.79E-04, 0.0, 0.0, 0.0, 0.0, 2.09E-04, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 1.39E-04};
	}

	/**
	 * The Destructor
	 */
	~A800HNeutronFluxHandler()
	{
	}
};
// end class A800HNeutronFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
