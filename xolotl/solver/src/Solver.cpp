// Includes
#include <xolotl/solver/Solver.h>

namespace xolotl
{
namespace solver
{
handler::ISolverHandler* Solver::solverHandler = nullptr;

Solver::Solver(handler::ISolverHandler& _solverHandler,
	std::shared_ptr<perf::IHandlerRegistry> registry) :
	optionsString(""), handlerRegistry(registry)
{
	solverHandler = &_solverHandler;
}

void
Solver::setCommandLineOptions(std::string arg)
{
	// Keep the arguments
	optionsString = arg;
}

} /* end namespace solver */
} /* end namespace xolotl */
