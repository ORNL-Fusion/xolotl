#ifndef INTERFACE_H
#define INTERFACE_H
#include <PetscSolver.h>

/**
 * Class defining the method to be coupled to another code through MOOSEApps
 */
class XolotlInterface {
private:
	/**
	 * The solver
	 */
	std::shared_ptr<xolotlSolver::PetscSolver> solver;

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
	void initializeXolotl(int argc, char **argv, MPI_Comm comm = MPI_COMM_WORLD,
			bool isStandalone = true);

	/**
	 * Set the final time and the dt.
	 *
	 * @param finalTime The wanted final time
	 * @param dt The wanted max time step
	 */
	void setTimes(double finalTime, double dt);

	/**
	 * Run the PETSc solve
	 */
	void solveXolotl();

	/**
	 * Get the local Xe rate that needs to be passed
	 *
	 * @return The local vector of rates
	 */
	std::vector<std::vector<std::vector<double> > > * getLocalXeRate();

	/**
	 * Get the local Xe conc
	 *
	 * @return The local vector of concs
	 */
	std::vector<std::vector<std::vector<double> > > * getLocalXeConc();

	/**
	 * Get the local Xe rate that needs to be passed
	 *
	 * @param xs, xm The start and width in the X direction on the local MPI process
	 * @param Mx The total width in the X direction
	 * @param ys, ym The start and width in the Y direction on the local MPI process
	 * @param My The total width in the Y direction
	 * @param zs, zm The start and width in the Z direction on the local MPI process
	 * @param Mz The total width in the Z direction
	 */
	void getLocalCoordinates(int &xs, int &xm, int &Mx, int &ys, int &ym,
			int &My, int &zs, int &zm, int &Mz);

	/**
	 * Set the location of one GB grid point.
	 *
	 * @param i, j, k The coordinate of the GB
	 */
	void setGBLocation(int i, int j = 0, int k = 0);

	/**
	 * Get the TS from the solver
	 *
	 * @return The TS
	 */
	TS & getTS() {
		return solver->getTS();
	}

	/**
	 * Get the grid information
	 *
	 * @param hy The spacing in the Y direction
	 * @param hz The spacing in the Z direction
	 * @return The grid in the X direction
	 */
	std::vector<double> getGridInfo(double &hy, double &hz);

	/**
	 * Finalize the solve
	 *
	 * @param isStandalone To know is Xolotl is used as a subcomponent of another code
	 */
	void finalizeXolotl(bool isStandalone = true);

};
// End class interface
#endif
