// Includes
#include <Solver.h>

using namespace xolotlCore;

namespace xolotlSolver {

ISolverHandler* Solver::solverHandler = nullptr;

Solver::Solver(ISolverHandler& _solverHandler,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		numCLIArgs(0), CLIArgs(NULL), handlerRegistry(registry) {

	solverHandler = &_solverHandler;
}

void Solver::setCommandLineOptions(std::string argv) {
	// Keep the arguments
	optionsString = argv;
}

} /* end namespace xolotlSolver */
