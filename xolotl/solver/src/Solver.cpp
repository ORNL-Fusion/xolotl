// Includes
#include <xolotl/core/network/INetworkHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>
#include <xolotl/factory/temperature/TemperatureHandlerFactory.h>
#include <xolotl/perf/PerfHandlerRegistry.h>
#include <xolotl/solver/Solver.h>

namespace xolotl
{
namespace solver
{
Solver::Solver(
	const options::IOptions& options, SolverHandlerGenerator handlerGenerator) :
	network(factory::network::NetworkHandlerFactory::get()
				.generate(options)
				->getNetwork()),
	materialHandler(
		factory::material::MaterialHandlerFactory::get().generate(options)),
	temperatureHandler(
		factory::temperature::TemperatureHandlerFactory::get().generate(
			options)),
	solverHandler(handlerGenerator(*network)),
	perfHandler(perf::PerfHandlerRegistry::get())

{
	assert(solverHandler);
	solverHandler->initializeHandlers(
		materialHandler, temperatureHandler, options);
}

Solver::Solver(const std::shared_ptr<handler::ISolverHandler>& _solverHandler,
	const std::shared_ptr<perf::IPerfHandler>& _perfHandler) :
	optionsString(""),
	solverHandler(_solverHandler),
	perfHandler(_perfHandler)
{
}

void
Solver::setCommandLineOptions(std::string arg)
{
	// Keep the arguments
	optionsString = arg;
}

} /* end namespace solver */
} /* end namespace xolotl */
