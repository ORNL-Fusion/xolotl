// Includes
#include <xolotl/core/network/INetworkHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>
#include <xolotl/factory/perf/PerfHandlerFactory.h>
#include <xolotl/factory/temperature/TemperatureHandlerFactory.h>
#include <xolotl/solver/Solver.h>

namespace xolotl
{
namespace solver
{
Solver::Solver(
	const options::IOptions& options, SolverHandlerGenerator handlerGenerator) :
	perfHandler(factory::perf::PerfHandlerFactory::get(perf::loadPerfHandlers)
					.generate(options)),
	initTimer([](auto&& timer) {
		timer->start();
		return timer;
	}(perfHandler->getTimer("Initialization"))),
	network(factory::network::NetworkHandlerFactory::get(
		core::network::loadNetworkHandlers)
				.generate(options)
				->getNetwork()),
	materialHandler(
		factory::material::MaterialHandlerFactory::get().generate(options)),
	temperatureHandler(
		factory::temperature::TemperatureHandlerFactory::get().generate(
			options)),
	solverHandler(handlerGenerator(*network, *perfHandler)),
    checkpointFile(options.getCheckpointFilePath())
{
	assert(solverHandler);
	solverHandler->initializeHandlers(
		materialHandler.get(), temperatureHandler.get(), options);
}

Solver::Solver(const std::shared_ptr<handler::ISolverHandler>& _solverHandler) :
	optionsString(""),
	solverHandler(_solverHandler),
	perfHandler(_solverHandler->getPerfHandler()),
    checkpointFile("")
{
}

void
Solver::setCommandLineOptions(std::string arg)
{
	// Keep the arguments
	optionsString = arg;
}

void
Solver::setExternalControlStep(std::size_t step)
{
    assert(monitor);
    monitor->setExternalControlStep(step);
}
} /* end namespace solver */
} /* end namespace xolotl */
