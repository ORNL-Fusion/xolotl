// Includes
#include <Solver.h>

using namespace xolotlCore;

namespace xolotlSolver {

// Allocate the static solver handler
ISolverHandler *Solver::solverHandler;

Solver::Solver(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
	numCLIArgs(0), CLIArgs(NULL), handlerRegistry(registry) {
}

void Solver::setCommandLineOptions(int argc, char **argv) {
	// Keep the arguments
	numCLIArgs = argc;
	CLIArgs = argv;
}

} /* end namespace xolotlSolver */
