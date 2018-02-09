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
bool initializeDimension(xolotlCore::Options &options,
		xolotlCore::IReactionNetwork& network);

/**
 * Access the created solver handler.
 *
 *  @return The solver handler.
 */
xolotlSolver::ISolverHandler& getSolverHandler();

}

#endif /* SOLVERHANDLERFACTORY_H */
