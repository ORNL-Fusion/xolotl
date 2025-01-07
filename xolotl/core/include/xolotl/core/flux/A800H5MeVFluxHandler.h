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
		this->fluxI = {0, 0.905167886, 0.001195838, 0.000597919, 0.000448439,
				5.97919E-05, 0.000104636, 0.000134532, 0.0, 0.0,
			3.28855E-05, 0, 0, 0.0, 0, 0, 0, 0, 0.0, 0, 0, 0, 0,
			0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0};

		this->highFluxI = {0, 0.028529737, 0.005994826, 0.002314166, 0.001583792,
			0.000565448, 0.000571995, 0.000505504, 0.000258115, 0.000182199,
			0.000242672, 0, 0, 0.000197382, 0, 0, 0, 0, 0.000212566, 0, 0, 0, 0,
			0.000197382, 0, 0, 0, 0, 0.000121466, 0, 0, 0, 0, 2.27749E-05, 0, 0, 0,
			0, 2.27749E-05, 0, 0, 0, 0, 2.27749E-05};

		this->fluxV = {0, 0.967460772162609, 0.00218701856655623, 0.000662817314779118, 0.000515291712103164,
				0.000347909447205322, 0.000183883054713605, 0.00035406144642229, 0.000105379595883434, 0.000111672181790613,
				0.000108889333299122, 0, 0, 0.00019507016312254, 0, 0, 0, 0, 4.71943943038404E-05, 0, 0, 0, 0,
				9.43887886076808E-05, 0, 0, 0, 0, 4.09018083966617E-05, 0, 0, 0, 0, 5.66332731646085E-05, 0, 0,
			0, 0, 2.83166365823043E-05, 0, 0, 0, 0, 1.88777577215362E-05, 0, 0, 0, 0, 4.40481013502511E-05,
			0, 0, 0, 0, 2.51703436287149E-05, 0, 0, 0, 0, 1.88777577215362E-05, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 1.25851718143574E-05};
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
