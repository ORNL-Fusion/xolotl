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
		this->fluxI = {0.0, 0.969022101, 0.003964995, 0.001828305, 0.001336963,
			0.000258214, 0.000357645, 0.000407942, 5.82507E-05, 4.11181E-05,
			0.000126598, 0.0, 0.0, 4.45446E-05, 0.0, 0.0, 0.0, 0.0, 4.79711E-05,
			0.0, 0.0, 0.0, 0.0, 4.45446E-05, 0.0, 0.0, 0.0, 0.0, 2.74121E-05,
			0.0, 0.0, 0.0, 0.0, 5.13977E-06, 0.0, 0.0, 0.0, 0.0, 5.13977E-06,
			0.0, 0.0, 0.0, 0.0, 5.13977E-06};

		this->fluxV = {0.0, 0.986715421, 0.001155166, 0.000397685, 0.000281691,
			0.000169486, 8.28486E-05, 0.000278385, 5.85161E-05, 5.99362E-05,
			7.75023E-05, 0.0, 0.0, 4.40228E-05, 0.0, 0.0, 0.0, 0.0, 1.06507E-05,
			0.0, 0.0, 0.0, 0.0, 2.13014E-05, 0.0, 0.0, 0.0, 0.0, 9.23059E-06,
			0.0, 0.0, 0.0, 0.0, 1.27808E-05, 0.0, 0.0, 0.0, 0.0, 6.39041E-06,
			0.0, 0.0, 0.0, 0.0, 4.26027E-06, 0.0, 0.0, 0.0, 0.0, 9.94064E-06,
			0.0, 0.0, 0.0, 0.0, 5.68036E-06, 0.0, 0.0, 0.0, 0.0, 4.26027E-06,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.84018E-06};
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
