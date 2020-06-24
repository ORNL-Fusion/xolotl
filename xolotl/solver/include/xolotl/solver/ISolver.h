#ifndef ISOLVER_H
#define ISOLVER_H

// Includes
#include <map>

#include <xolotl/options/Options.h>
#include <xolotl/solver/handler/ISolverHandler.h>

namespace xolotl
{
namespace solver
{
/**
 * Realizations of this interface are responsible for solving the set of
 * advection-diffusion-reaction equations that govern the evolution of
 * clusters in plasma surface interactions.
 */
class ISolver
{
public:
	/**
	 * The destructor
	 */
	virtual ~ISolver()
	{
	}

	/**
	 * This operation transfers the input arguments passed to the program on
	 * startup to the solver. These options are static options specified at
	 * the start of the program whereas the options passed to setOptions() may
	 * change.
	 *
	 * @param arg The string containing all the solver options
	 */
	virtual void
	setCommandLineOptions(std::string arg) = 0;

	/**
	 * This operation sets the run-time options of the solver. The map is a set
	 * of key-value std::string pairs that are interpreted by the solver. These
	 * options may change during execution, but it is up to Solvers to monitor
	 * the map for changes and they may do so at their discretion.
	 * @param options The set of options as key-value pairs with option names
	 * for keys and associated values mapped to those keys. A relevant example
	 * is "startTime" and "0.01" where both are of type std::string.
	 */
	virtual void
	setOptions(const std::map<std::string, std::string>& options) = 0;

	/**
	 * This operation sets up the mesh that will be used by the solver and
	 * initializes the data on that mesh. This operation will throw an exception
	 * of type std::string if the mesh can not be setup.
	 */
	virtual void
	setupMesh() = 0;

	/**
	 * This operation performs all necessary initialization for the solver
	 * possibly including but not limited to setting up MPI and loading initial
	 * conditions. If the solver can not be initialized, this operation will
	 * throw an exception of type std::string.
	 *
	 * @param isStandalone To know is Xolotl is used as a subcomponent of
	 * another code
	 */
	virtual void
	initialize() = 0;

	/**
	 * This operation sets the wanted final time and max time step allowed.
	 *
	 * @param finalTime The wanted final time.
	 * @param dt The wanted max time step.
	 */
	virtual void
	setTimes(double finalTime, double dt) = 0;

	/**
	 * This operation resets the concentrations where the GB are located.
	 */
	virtual void
	initGBLocation() = 0;

	/**
	 * This operation get the concentration vector with the ids.
	 *
	 * @return The concentration vector
	 */
	virtual std::vector<
		std::vector<std::vector<std::vector<std::pair<int, double>>>>>
	getConcVector() = 0;

	/**
	 * This operation sets the concentration vector in the current state of the
	 * simulation.
	 *
	 * @param The concentration vector
	 */
	virtual void
	setConcVector(std::vector<std::vector<
			std::vector<std::vector<std::pair<int, double>>>>>& concVector) = 0;

	/**
	 * Get the current dt.
	 *
	 * @return  The current time step
	 */
	virtual double
	getCurrentDt() = 0;

	/**
	 * Set the current time and dt.
	 *
	 * @param currentTime The time
	 * @param currentDt The current time step
	 */
	virtual void
	setCurrentTimes(double currentTime, double currentDt) = 0;

	/**
	 * This operation directs the Solver to perform the solve. If the solve
	 * fails, it will throw an exception of type std::string.
	 */
	virtual void
	solve() = 0;

	/**
	 * This operation checks the solver convergence status to decide whether
	 * it converged or not.
	 *
	 * @return true if it converged
	 */
	virtual bool
	getConvergenceStatus() = 0;

	/**
	 * This operation performs all necessary finalization for the solver
	 * including but not limited to cleaning up memory, finalizing MPI and
	 * printing diagnostic information. If the solver can not be finalized,
	 * this operation will throw an exception of type std::string.
	 */
	virtual void
	finalize() = 0;

	/**
	 * Get the TS.
	 *
	 * @return The TS
	 */
	virtual TS&
	getTS() = 0;

	/**
	 * This operation returns the most recent time that solver converged.
	 */
	virtual double
	getXolotlTime() = 0;
};
// end class ISolver

} /* namespace solver */
} /* namespace xolotl */
#endif
