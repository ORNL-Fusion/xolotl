#ifndef SOLVERHANDLERFACTORY_H
#define SOLVERHANDLERFACTORY_H

#include <memory>

#include <xolotl/options/Options.h>
#include <xolotl/solver/handler/ISolverHandler.h>

namespace xolotl
{
namespace factory
{
namespace solver
{
/**
 * Build the desired type of solver.
 *
 * @param opts Options for the program
 * @return True if the solver handler was created successfully.
 */
bool
initializeDimension(
	const options::Options& opts, core::network::IReactionNetwork& network);

/**
 * Access the created solver handler.
 *
 *  @return The solver handler.
 */
::xolotl::solver::handler::ISolverHandler&
getSolverHandler();

/**
 * Destroy the solver handler
 */
void
destroySolverHandler();

} // namespace solver
} // namespace factory
} // namespace xolotl

#endif /* SOLVERHANDLERFACTORY_H */
