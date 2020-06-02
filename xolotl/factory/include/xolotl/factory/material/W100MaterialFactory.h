#ifndef W100MATERIALHANDLERFACTORY_H
#define W100MATERIALHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/material/MaterialFactory.h>
#include <xolotl/core/flux/W100FitFluxHandler.h>
#include <xolotl/core/advection/W100AdvectionHandler.h>
#include <xolotl/core/modified/W100TrapMutationHandler.h>

namespace xolotl {
namespace factory {
namespace material {

/**
 * Subclass of MaterialFactory for a (100) oriented tungsten material.
 */
class W100MaterialFactory: public MaterialFactory {
public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W100MaterialFactory(const options::Options &opts) :
			MaterialFactory(opts) {
		theFluxHandler = std::make_shared<core::flux::W100FitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<core::advection::W100AdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				core::modified::W100TrapMutationHandler>();
		theNucleationHandler = std::make_shared<
				core::modified::DummyNucleationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~W100MaterialFactory() {
	}
};

} // end namespace material
} // end namespace factory
} // end namespace xolotl

#endif // W100MATERIALHANDLERFACTORY_H
