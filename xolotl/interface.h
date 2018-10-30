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
	 */
	std::shared_ptr<xolotlSolver::PetscSolver> initializeXolotl(int argc,
			char **argv, MPI_Comm comm = MPI_COMM_WORLD);

	/**
	 * Run the PETSc solve
	 */
	void solveXolotl(std::shared_ptr<xolotlSolver::PetscSolver> solver);

	/**
	 * Get the retention vector copy at the end of solve
	 *
	 * @return The vector
	 */
	virtual std::vector<double> getCopyRetention(std::shared_ptr<xolotlSolver::PetscSolver> solver);

	/**
	 * Get the retention vector pointer at the end of solve
	 *
	 * @return The vector
	 */
	virtual std::vector<double>* getPointerRetention(std::shared_ptr<xolotlSolver::PetscSolver> solver);

	/**
	 * Print the retention vector copy at the end of solve
	 */
	virtual void printRetention(std::shared_ptr<xolotlSolver::PetscSolver> solver);

	/**
	 * Finalize the solve
	 */
	void finalizeXolotl(std::shared_ptr<xolotlSolver::PetscSolver> solver);

};
// End class interface
#endif
