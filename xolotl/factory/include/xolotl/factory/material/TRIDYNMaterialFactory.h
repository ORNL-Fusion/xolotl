#ifndef TRIDYNMATERIALHANDLERFACTORY_H
#define TRIDYNMATERIALHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/material/MaterialFactory.h>
#include <xolotl/core/flux/TRIDYNFitFluxHandler.h>
#include <xolotl/core/advection/W100AdvectionHandler.h>
#include <xolotl/core/modifiedreaction/trapmutation/W100TrapMutationHandler.h>

namespace xolotlFactory {

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
	TRIDYNMaterialFactory(const xolotlCore::Options &options) :
			MaterialFactory(options) {
		theFluxHandler = std::make_shared<xolotlCore::TRIDYNFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::W100AdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::W100TrapMutationHandler>();
		theNucleationHandler = std::make_shared<
				xolotlCore::DummyNucleationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~TRIDYNMaterialFactory() {
	}
};

} // end namespace xolotlFactory

#endif // TRIDYNMATERIALHANDLERFACTORY_H
