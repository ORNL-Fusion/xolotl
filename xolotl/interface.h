#ifndef INTERFACE_H
#define INTERFACE_H
#include <PetscSolver.h>

/**
 * Class defining the method to be coupled to another code through MOOSEApps
 */
class XolotlInterface {

public:

	/**
	 * The default constructor
	 */
	XolotlInterface() {
	}

	/**
	 * The destructor
	 */
	~XolotlInterface() {
	}

	/**
	 * Print something
	 */
	void printSomething();

	/**
	 * Initialize all the options and handlers
	 *
	 * @param argc, argv The command line arguments
	 * @param MPI_Comm The communicator to use
	 * @param isStandalone To know is Xolotl is used as a subcomponent of another code
	 * @return The pointer to the solver
	 */
	std::shared_ptr<xolotlSolver::PetscSolver> initializeXolotl(int argc,
			char **argv, MPI_Comm comm = MPI_COMM_WORLD, bool isStandalone =
					true);

	/**
	 * Set the final time and the dt.
	 *
	 * @param solver The pointer to the solver
	 * @param finalTime The wanted final time
	 * @param dt The wanted max time step
	 */
	void setTimes(std::shared_ptr<xolotlSolver::PetscSolver> solver,
			double finalTime, double dt);

	/**
	 * Run the PETSc solve
	 *
	 * @param solver The pointer to the solver
	 */
	void solveXolotl(std::shared_ptr<xolotlSolver::PetscSolver> solver);

	/**
	 * Get the local Xe rate that needs to be passed
	 *
	 * @param solver The pointer to the solver
	 * @return The local vector of rates
	 */
	std::vector<std::vector<std::vector<double> > > * getLocalXeRate(
			std::shared_ptr<xolotlSolver::PetscSolver> solver);

	/**
	 * Get the local Xe rate that needs to be passed
	 *
	 * @param solver The pointer to the solver
	 * @param xs, xm The start and width in the X direction on the local MPI process
	 * @param Mx The total width in the X direction
	 * @param ys, ym The start and width in the Y direction on the local MPI process
	 * @param My The total width in the Y direction
	 * @param zs, zm The start and width in the Z direction on the local MPI process
	 * @param Mz The total width in the Z direction
	 */
	void getLocalCoordinates(std::shared_ptr<xolotlSolver::PetscSolver> solver,
			int &xs, int &xm, int &Mx, int &ys, int &ym, int &My, int &zs,
			int &zm, int &Mz);

	/**
	 * Finalize the solve
	 *
	 * @param solver The pointer to the solver
	 * @param isStandalone To know is Xolotl is used as a subcomponent of another code
	 */
	void finalizeXolotl(std::shared_ptr<xolotlSolver::PetscSolver> solver,
			bool isStandalone = true);

};
// End class interface
#endif
