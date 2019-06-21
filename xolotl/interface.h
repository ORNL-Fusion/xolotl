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
	 * @return The vector
	 */
	std::vector<std::vector<std::vector<double> > > getLocalXeRate();

	/**
	 * Set the local Xe rate
	 *
	 * @param The vector
	 */
	void setLocalXeRate(std::vector<std::vector<std::vector<double> > > rateVector);

	/**
	 * Get the previous Xe flux that needs to be passed
	 *
	 * @return The vector
	 */
	std::vector<std::vector<std::vector<double> > > getPreviousXeFlux();

	/**
	 * Set the previous Xe flux
	 *
	 * @param The vector
	 */
	void setPreviousXeFlux(std::vector<std::vector<std::vector<double> > > fluxVector);

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
	 * Reset the GB vector.
	 */
	void resetGBVector();

	/**
	 * Set the concentrations to 0.0 where the GBs are.
	 */
	void initGBLocation();

	/**
	 * Get the concentrations and their ids.
	 *
	 * @return The concentration vector from the current state of the simulation
	 */
	std::vector<std::vector<std::vector<std::vector<std::pair<int, double> > > > > getConcVector();

	/**
	 * Set the concentrations and their ids.
	 *
	 * @ param concVector A given state of the concentrations
	 */
	void setConcVector(std::vector<std::vector<std::vector<std::vector<std::pair<int, double> > > > > concVector);

	/**
	 * Get the previous time.
	 *
	 * @return The previous time
	 */
	double getPreviousTime();

	/**
	 * Set the previous time.
	 *
	 * @param time The previous time
	 */
	void setPreviousTime(double time);

	/**
	 * Get the current dt.
	 *
	 * @return  The current time step
	 */
	double getCurrentDt();

	/**
	 * Set the current time and dt.
	 *
	 * @param currentTime The time
	 * @param currentDt The current time step
	 */
	void setCurrentTimes(double currentTime, double currentDt);

	/**
	 * Get the number of Xe that went to the GB.
	 *
	 * @return The number of Xenon
	 */
	double getNXeGB();

	/**
	 * Set the number of Xe that went to the GB.
	 *
	 * @param nXe The number of Xenon
	 */
	void setNXeGB(double nXe);

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
