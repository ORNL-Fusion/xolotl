#include <xolotl/factory/reaction/IReactionHandlerFactory.h>
#include <xolotl/factory/reaction/PSIReactionHandlerFactory.h>
#include <xolotl/factory/reaction/NEReactionHandlerFactory.h>
#include <xolotl/factory/reaction/AlloyReactionHandlerFactory.h>
#include <xolotl/factory/reaction/FeReactionHandlerFactory.h>

namespace xolotlFactory {

std::shared_ptr<IReactionHandlerFactory> theReactionFactory;

std::shared_ptr<IReactionHandlerFactory> IReactionHandlerFactory::createNetworkFactory(
		const std::string& problemType) {
	// PSI case
	if (problemType == "W100" || problemType == "W110" || problemType == "W111"
			|| problemType == "W211" || problemType == "TRIDYN" || problemType == "Pulsed")
		theReactionFactory = std::make_shared<PSIReactionHandlerFactory>();
	// NE case
	else if (problemType == "Fuel")
		theReactionFactory = std::make_shared<NEReactionHandlerFactory>();
	// Alloy case
	else if (problemType == "800H")
		theReactionFactory = std::make_shared<AlloyReactionHandlerFactory>();
	// Fe case
	else if (problemType == "Fe")
		theReactionFactory = std::make_shared<FeReactionHandlerFactory>();
	// The type is not supported
	else {
		throw std::string(
				"\nThe problem type is not known: \"" + problemType + "\"");
	}

	return theReactionFactory;
}

void IReactionHandlerFactory::resetNetworkFactory()
{
    theReactionFactory.reset();
}

} // end namespace xolotlFactory
