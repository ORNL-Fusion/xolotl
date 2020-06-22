#include <cassert>
#include <fstream>
#include <iostream>

#include <mpi.h>

#include <xolotl/factory/solver/SolverHandlerFactory.h>
#include <xolotl/solver/handler/PetscSolver0DHandler.h>
#include <xolotl/solver/handler/PetscSolver1DHandler.h>
#include <xolotl/solver/handler/PetscSolver2DHandler.h>
#include <xolotl/solver/handler/PetscSolver3DHandler.h>

using namespace xolotl::solver::handler;

namespace xolotl
{
namespace factory
{
namespace solver
{
std::unique_ptr<ISolverHandler> theSolverHandler;

// Create the desired type of handler registry.
bool
initializeDimension(
	const options::Options& opts, core::network::IReactionNetwork& network)
{
	bool ret = true;

	// Get the wanted dimension
	int dim = opts.getDimensionNumber();

	// Switch on the dimension
	// TODO Once we have widespread C++14 support, use std::make_unique
	// instead of this two-step construction.
	ISolverHandler* rawSolverHandler = nullptr;
	switch (dim) {
	case 0:
		rawSolverHandler = new PetscSolver0DHandler(network);
		break;
	case 1:
		rawSolverHandler = new PetscSolver1DHandler(network);
		break;
	case 2:
		rawSolverHandler = new PetscSolver2DHandler(network);
		break;
	case 3:
		rawSolverHandler = new PetscSolver3DHandler(network);
		break;
	default:
		// The asked dimension is not good (e.g. -1, 4)
		throw std::string(
			"\nxolotlFactory: Bad dimension for the solver handler.");
	}
	assert(rawSolverHandler != nullptr);
	theSolverHandler = std::unique_ptr<ISolverHandler>(rawSolverHandler);

	return ret;
}

// Provide access to our handler registry.
ISolverHandler&
getSolverHandler()
{
	if (!theSolverHandler) {
		// We have not yet been initialized.
		throw std::string("\nxolotlFactory: solver requested but "
						  "it has not been initialized.");
	}

	return *theSolverHandler;
}

void
destroySolverHandler()
{
	if (theSolverHandler) {
		theSolverHandler.reset();
	}

	return;
}

} // end namespace solver
} // end namespace factory
} // end namespace xolotl
