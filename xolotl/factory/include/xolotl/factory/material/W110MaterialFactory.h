#ifndef W110MATERIALHANDLERFACTORY_H
#define W110MATERIALHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/material/MaterialFactory.h>
#include <xolotl/core/flux/W110FitFluxHandler.h>
#include <xolotl/core/advection/W110AdvectionHandler.h>
#include <xolotl/core/modified/W110TrapMutationHandler.h>

namespace xolotl {
namespace factory {
namespace material {

/**
 * Subclass of MaterialFactory for a (110) oriented tungsten material.
 */
class W110MaterialFactory: public MaterialFactory {

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W110MaterialFactory(const options::Options &opts) :
			MaterialFactory(opts) {
		theFluxHandler = std::make_shared<core::flux::W110FitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<core::advection::W110AdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				core::modified::W110TrapMutationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~W110MaterialFactory() {
	}
};

} // end namespace material
} // end namespace factory
} // end namespace xolotl

#endif // W110MATERIALHANDLERFACTORY_H
