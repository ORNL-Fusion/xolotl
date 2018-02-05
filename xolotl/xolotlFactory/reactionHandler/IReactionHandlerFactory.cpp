#include "IReactionHandlerFactory.h"
#include "PSIReactionHandlerFactory.h"
#include "NEReactionHandlerFactory.h"
#include "FeReactionHandlerFactory.h"

namespace xolotlFactory {

static std::shared_ptr<IReactionHandlerFactory> theReactionFactory;

std::shared_ptr<IReactionHandlerFactory> IReactionHandlerFactory::createNetworkFactory(
		const std::string& problemType) {
	// PSI case
	if (problemType == "W100" || problemType == "W110" || problemType == "W111"
			|| problemType == "W211" || problemType == "TRIDYN")
		theReactionFactory = std::make_shared<PSIReactionHandlerFactory>();
	// NE case
	else if (problemType == "Fuel")
		theReactionFactory = std::make_shared<NEReactionHandlerFactory>();
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

} // end namespace xolotlFactory
