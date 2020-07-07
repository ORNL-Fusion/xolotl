#ifndef SOLVER_H
#define SOLVER_H

// Includes
#include <xolotl/solver/ISolver.h>

namespace xolotl
{
namespace solver
{
/**
 * This class and its subclasses realize the ISolver interface to solve the
 * advection-diffusion-reaction problem with currently supported solvers.
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

	//! Static global reference to current solver handler
	static handler::ISolverHandler* staticSolverHandler;

public:
	using SolverHandlerGenerator =
		std::function<std::shared_ptr<handler::ISolverHandler>(
			core::network::IReactionNetwork&)>;

	/**
	 * Default constructor, deleted because we must have arguments to construct.
	 */
	Solver() = delete;

	Solver(const options::Options& options,
		SolverHandlerGenerator handlerGenerator);

	//! Constuct a solver.
	Solver(handler::ISolverHandler& _solverHandler,
		std::shared_ptr<perf::IHandlerRegistry> registry);

	//! The Destructor
	virtual ~Solver(){};

	/**
	 * This operation transfers the input arguments passed to the program on
	 * startup to the solver. These options are static options specified at
	 * the start of the program whereas the options passed to setOptions() may
	 * change.
	 * @param arg The string of arguments
	 */
	void
	setCommandLineOptions(std::string arg);

	/**
	 * This operation returns the solver handler for this solver. This
	 * operation is only for use by solver code and is not part of the
	 * ISolver interface.
	 * @return The advection handler for this solver
	 */
	static handler::ISolverHandler&
	getSolverHandler()
	{
		return *staticSolverHandler;
	}

protected:
	/**
	 * The performance handler registry that will be used
	 * for this class.
	 */
	std::shared_ptr<perf::IHandlerRegistry> handlerRegistry;
};
// end class Solver

} /* namespace solver */
} /* namespace xolotl */
#endif
