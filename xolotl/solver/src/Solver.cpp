// Includes
#include <xolotl/core/network/INetworkHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>
#include <xolotl/factory/temperature/TemperatureHandlerFactory.h>
#include <xolotl/solver/Solver.h>

namespace xolotl
{
namespace solver
{
handler::ISolverHandler* Solver::staticSolverHandler = nullptr;

Solver::Solver(
	const options::Options& options, SolverHandlerGenerator handlerGenerator) :
	network(factory::network::NetworkHandlerFactory::get()
				.generateNetworkHandler(options)
				->getNetwork()),
	materialHandler(factory::material::MaterialHandlerFactory::get()
						.generateMaterialHandler(options)),
	temperatureHandler(factory::temperature::TemperatureHandlerFactory::get()
						   .generateTemperatureHandler(options)),
	solverHandler(handlerGenerator(*network)),
	handlerRegistry(perf::getHandlerRegistry())

{
	assert(solverHandler);
	solverHandler->initializeHandlers(
		materialHandler, temperatureHandler, options);
	staticSolverHandler = solverHandler.get();
}

Solver::Solver(handler::ISolverHandler& _solverHandler,
	std::shared_ptr<perf::IHandlerRegistry> registry) :
	optionsString(""), handlerRegistry(registry)
{
	staticSolverHandler = &_solverHandler;
}

void
Solver::setCommandLineOptions(std::string arg)
{
	// Keep the arguments
	optionsString = arg;
}

} /* end namespace solver */
} /* end namespace xolotl */
