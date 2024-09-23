#ifndef SOLVER_H
#define SOLVER_H

// Includes
#include <xolotl/options/IOptions.h>
#include <xolotl/perf/IPerfHandler.h>
#include <xolotl/perf/ITimer.h>
#include <xolotl/solver/ISolver.h>
#include <xolotl/solver/monitor/IMonitor.h>

namespace xolotl
{
namespace solver
{
/**
 * This class and its subclasses realize the ISolver interface to solve the
 * diffusion-reaction problem with currently supported solvers.
 */
class Solver : public ISolver
{
protected:
	//! The string of option
	std::string optionsString;

	//! The perf handler
	std::shared_ptr<perf::IPerfHandler> perfHandler;

	//! The initialization timer
	std::shared_ptr<perf::ITimer> initTimer;

	//! The network
	std::shared_ptr<core::network::IReactionNetwork> network;

	//! The material handler
	std::shared_ptr<core::material::IMaterialHandler> materialHandler;

	//! The temperature handler
	std::shared_ptr<core::temperature::ITemperatureHandler> temperatureHandler;

	//! The original solver handler.
	std::shared_ptr<handler::ISolverHandler> solverHandler;

	//! The checkpoint file name
	std::string checkpointFile;

	//! The monitor
	std::shared_ptr<monitor::IMonitor> monitor;

public:
	using SolverHandlerGenerator =
		std::function<std::shared_ptr<handler::ISolverHandler>(
			core::network::IReactionNetwork&, perf::IPerfHandler&)>;

	/**
	 * Default constructor, deleted because we must have arguments to construct.
	 */
	Solver() = delete;

	Solver(const options::IOptions& options,
		SolverHandlerGenerator handlerGenerator);

	//! Constuct a solver.
	Solver(const std::shared_ptr<handler::ISolverHandler>& _solverHandler);

	//! The Destructor
	virtual ~Solver()
	{
	}

	/**
	 * \see ISolver.h
	 */
	void
	setCommandLineOptions(std::string arg);

	/**
	 * \see ISolver.h
	 */
	void
	setExternalControlStep(std::size_t step) override;

	/**
	 * @return The solver handler for this solver
	 */
	handler::ISolverHandler*
	getSolverHandler()
	{
		return solverHandler.get();
	}
};
// end class Solver

} /* namespace solver */
} /* namespace xolotl */
#endif
