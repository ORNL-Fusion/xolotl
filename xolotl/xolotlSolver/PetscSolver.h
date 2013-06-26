#ifndef PETSCSOLVER_H
#define PETSCSOLVER_H

// Includes
#include "ISolver.h"
#include <petscsys.h>
#include <petscdmda.h>
#include <petscts.h>

namespace xolotlSolver {

/**
 * This class realizes the ISolver interface to solve the
 * advection-diffusion-reaction problem with the Petsc solvers from Argonne
 * National Laboratory.
 */
class PetscSolver : public ISolver {

public:

	//! The Constructor
	PetscSolver();

	//! The Destructor
	~PetscSolver();

	/**
	 * This operation transfers the input arguments passed to the program on
	 * startup to the solver. These options are static options specified at
	 * the start of the program whereas the options passed to setOptions() may
	 * change.
	 * @param argc The number of command line arguments
	 * @param argv The array of command line arguments
	 */
	void setCommandLineOptions(int argc, char **argv);
	
	/**
	 * This operation sets the PSIClusterNetworkLoader that should be used by
	 * the ISolver to load the ReactionNetwork.
	 * @param networkLoader The PSIClusterNetworkLoader that will load the
	 * network.
	 */
	void setNetworkLoader(std::shared_ptr<PSIClusterNetworkLoader> networkLoader);

	/**
	 * This operation sets the run-time options of the solver. The map is a set
	 * of key-value std::string pairs that are interpreted by the solver. These
	 * options may change during execution, but it is up to Solvers to monitor
	 * the map for changes and they may do so at their discretion.
	 * @param options The set of options as key-value pairs with option names
	 * for keys and associated values mapped to those keys. A relevant example
	 * is "startTime" and "0.01" where both are of type std::string.
	 */
	void setOptions(std::map<std::string, std::string> options);

	/**
	 * This operation sets up the mesh that will be used by the solver and
	 * initializes the data on that mesh. This operation will throw an exception
	 * of type std::string if the mesh can not be setup.
	 */
	void setupMesh();

	/**
	 * This operation performs all necessary initialization for the solver
	 * possibly including but not limited to setting up MPI and loading initial
	 * conditions. If the solver can not be initialized, this operation will
	 * throw an exception of type std::string.
	 */
	void initialize();

	/**
	 * This operation directs the Solver to perform the solve. If the solve
	 * fails, it will throw an exception of type std::string.
	 */
	void solve();

	/**
	 * This operation performs all necessary finalization for the solver
	 * including but not limited to cleaning up memory, finalizing MPI and
	 * printing diagnostic information. If the solver can not be finalized,
	 * this operation will throw an exception of type std::string.
	 */
	void finalize();
	
private:
	//! The number command line arguments
	int numCLIArgs;

	//! The command line arguments
	char **CLIArgs;
	
	std::shared_ptr<PSIClusterNetworkLoader> networkLoader;
	
}; //end class PetscSolver

} /* end namespace xolotlSolver */
#endif
