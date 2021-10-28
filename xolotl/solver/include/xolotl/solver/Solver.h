#ifndef SOLVER_H
#define SOLVER_H

// Includes
#include <xolotl/options/IOptions.h>
#include <xolotl/perf/IPerfHandler.h>
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

	//! The network
	std::shared_ptr<core::network::IReactionNetwork> network;

	//! The material handler
	std::shared_ptr<core::material::IMaterialHandler> materialHandler;

	//! The temperature handler
	std::shared_ptr<core::temperature::ITemperatureHandler> temperatureHandler;

	//! The original solver handler.
	std::shared_ptr<handler::ISolverHandler> solverHandler;

	//! The monitor
	std::shared_ptr<monitor::IMonitor> monitor;

public:
	using SolverHandlerGenerator =
		std::function<std::shared_ptr<handler::ISolverHandler>(
			core::network::IReactionNetwork&)>;

	/**
	 * Default constructor, deleted because we must have arguments to construct.
	 */
	Solver() = delete;

	Solver(const options::IOptions& options,
		SolverHandlerGenerator handlerGenerator);

	//! Constuct a solver.
	Solver(const std::shared_ptr<handler::ISolverHandler>& _solverHandler,
		const std::shared_ptr<perf::IPerfHandler>& _perfHandler);

	//! The Destructor
	virtual ~Solver(){};

	/**
	 * \see ISolver.h
	 */
	void
	setCommandLineOptions(std::string arg);

	/**
	 * @return The solver handler for this solver
	 */
	std::shared_ptr<handler::ISolverHandler>
	getSolverHandler()
	{
		return solverHandler;
	}

protected:
	/**
	 * The performance handler registry that will be used
	 * for this class.
	 */
	std::shared_ptr<perf::IPerfHandler> perfHandler;
};
// end class Solver

} /* namespace solver */
} /* namespace xolotl */
#endif
