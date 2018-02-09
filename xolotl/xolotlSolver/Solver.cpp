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

void Solver::setCommandLineOptions(int argc, char **argv) {
	// Keep the arguments
	numCLIArgs = argc;
	CLIArgs = argv;
}

} /* end namespace xolotlSolver */
