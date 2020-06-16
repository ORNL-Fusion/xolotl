#ifndef TRIDYNMATERIALHANDLERFACTORY_H
#define TRIDYNMATERIALHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/material/MaterialFactory.h>
#include <xolotl/core/flux/TRIDYNFitFluxHandler.h>
#include <xolotl/core/advection/W100AdvectionHandler.h>
#include <xolotl/core/modified/W100TrapMutationHandler.h>

namespace xolotl {
namespace factory {
namespace material {

/**
 * Subclass of MaterialFactory for a (100) oriented tungsten material with
 * a TRIDYN input file.
 */
class TRIDYNMaterialFactory: public MaterialFactory {
public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	TRIDYNMaterialFactory(const options::Options &opts) :
			MaterialFactory(opts) {
		theFluxHandler = std::make_shared<core::flux::TRIDYNFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<core::advection::W100AdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				core::modified::W100TrapMutationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~TRIDYNMaterialFactory() {
	}
};

} // end namespace material
} // end namespace factory
} // end namespace xolotl

#endif // TRIDYNMATERIALHANDLERFACTORY_H
