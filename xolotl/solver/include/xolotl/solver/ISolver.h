#ifndef ISOLVER_H
#define ISOLVER_H

// Includes
#include <map>

#include <xolotl/solver/handler/ISolverHandler.h>

namespace xolotl
{
namespace solver
{
/**
 * Realizations of this interface are responsible for solving the set of
 * diffusion-reaction equations that govern the evolution of
 * cluster concentrations.
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
	 * startup to the solver.
	 *
	 * @param arg The string containing all the solver options
	 */
	virtual void
	setCommandLineOptions(std::string arg) = 0;

	/**
	 * This operation performs all necessary initialization for the solver.
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
	 * fails, it will throw an exception of type std::runtime_error.
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
	 * This operation performs all necessary finalization for the solver.
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
	 *
	 * @return The time
	 */
	virtual double
	getXolotlTime() = 0;
};
// end class ISolver

} /* namespace solver */
} /* namespace xolotl */
#endif
