// Includes
#include <Solver.h>

using namespace xolotlCore;

namespace xolotlSolver {

ISolverHandler *Solver::solverHandler = nullptr;

Solver::Solver(ISolverHandler &_solverHandler,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		optionsString(""), handlerRegistry(registry) {

	solverHandler = &_solverHandler;
}

void Solver::setCommandLineOptions(std::string arg) {
	// Keep the arguments
	optionsString = arg;
}

} /* end namespace xolotlSolver */
