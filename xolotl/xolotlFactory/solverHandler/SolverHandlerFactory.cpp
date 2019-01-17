#include <cassert>
#include <SolverHandlerFactory.h>
#include <PetscSolver0DHandler.h>
#include <PetscSolver1DHandler.h>
#include <PetscSolver2DHandler.h>
#include <PetscSolver3DHandler.h>
#include <fstream>
#include <iostream>
#include <mpi.h>

namespace xolotlFactory {

static std::unique_ptr<xolotlSolver::ISolverHandler> theSolverHandler;

// Create the desired type of handler registry.
bool initializeDimension(const xolotlCore::Options &options,
		xolotlCore::IReactionNetwork& network) {

	bool ret = true;

	// Get the wanted dimension
	int dim = options.getDimensionNumber();

	// Switch on the dimension
	// TODO Once we have widespread C++14 support, use std::make_unique
	// instead of this two-step construction.
	xolotlSolver::ISolverHandler* rawSolverHandler = nullptr;
	switch (dim) {
	case 0:
		rawSolverHandler = new xolotlSolver::PetscSolver0DHandler(network);
		break;
	case 1:
		rawSolverHandler = new xolotlSolver::PetscSolver1DHandler(network);
		break;
	case 2:
		rawSolverHandler = new xolotlSolver::PetscSolver2DHandler(network);
		break;
	case 3:
		rawSolverHandler = new xolotlSolver::PetscSolver3DHandler(network);
		break;
	default:
		// The asked dimension is not good (e.g. -1, 4)
		throw std::string(
				"\nxolotlFactory: Bad dimension for the solver handler.");
	}
	assert(rawSolverHandler != nullptr);
	theSolverHandler = std::unique_ptr<xolotlSolver::ISolverHandler>(
			rawSolverHandler);

	// Take care of the free GB option for NE
	if (options.getMaterial() == "Fuel")
		rawSolverHandler->setGBFileName(options.getGbString());

	return ret;
}

// Provide access to our handler registry.
xolotlSolver::ISolverHandler& getSolverHandler() {
	if (!theSolverHandler) {
		// We have not yet been initialized.
		throw std::string("\nxolotlFactory: solver requested but "
				"it has not been initialized.");
	}

	return *theSolverHandler;
}

} // end namespace xolotlFactory

