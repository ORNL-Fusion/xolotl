#ifndef SOLVERHANDLERFACTORY_H
#define SOLVERHANDLERFACTORY_H

#include <memory>
#include <ISolverHandler.h>
#include <Options.h>

namespace xolotlFactory {

/**
 * Build the desired type of solver.
 *
 * @param options Options for the program
 * @return True if the solver handler was created successfully.
 */
bool initializeDimension(const xolotlCore::Options &options,
		xolotlCore::experimental::IReactionNetwork& network);

/**
 * Access the created solver handler.
 *
 *  @return The solver handler.
 */
std::unique_ptr<xolotlSolver::ISolverHandler> getSolverHandler();

/**
 * Destroy the solver handler
 */
void destroySolverHandler();

}

#endif /* SOLVERHANDLERFACTORY_H */
