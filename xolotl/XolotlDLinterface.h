#ifndef XOLOTLDLINTERFACE_H
#define XOLOTLDLINTERFACE_H

#include <PetscSolver.h>
// #include <tuple>
//#include <MPIUtils.h>

/**
 * Class defining the method to be coupled to another code through MOOSEApps
 */
class XolotlDLinterface {

public:

	/**
	 * The default constructor
	 */
	//XolotlInterface() {
	//}
	XolotlDLinterface(){
	}

	/**
	 * The destructor
	 */
	//~XolotlInterface() {
	//}
	virtual ~XolotlDLinterface(){};

	/**
	 * Print something
	 */
	virtual void printSomething(){};

	/**
	 * Initialize all the options and handlers
	 */
	virtual std::shared_ptr<xolotlSolver::PetscSolver> initializeXolotl(int argc, char **argv, MPI_Comm comm = MPI_COMM_WORLD, bool isStandalone = true){
		std::shared_ptr<xolotlSolver::PetscSolver> a;
		return a;
	};

	virtual void getXolotlGlobalGridInfo(int *dim, bool *regulargrid, int *nx, int *ny, int *nz, double *hx, double *hy, double *hz,
				double *lx, double *ly, double *lz,
				int argc, char **argv){};

	virtual std::vector<double> getXolotlXgrid(std::shared_ptr<xolotlSolver::PetscSolver> solver){
		std::vector<double> a = {0.0};
		return a;
	};

	/**
 	 * Set the final time and the dt.
 	 *
	 * @param finalTime The wanted final time.
 	 * @param dt The wanted max time step.
 	 */
	virtual void setTimes(std::shared_ptr<xolotlSolver::PetscSolver> solver,
			double finalTime, double dt){};

	/**
 	 * Get the most recent Xolotl time that the solver converged
 	 */
	virtual double getXolotlTimeInterface(std::shared_ptr<xolotlSolver::PetscSolver> solver){
		return 0;
	};

	/**
	 * Run the PETSc solve
	 */
	virtual void solveXolotl(std::shared_ptr<xolotlSolver::PetscSolver> solver){};

	/**
	 * Get the localNE that needs to be passed
	 *
	 * @param solver The pointer to the solver
	 * @return The local vector of rates
	 * localNE[0][i][j][k] : XeRate
	 * localNE[1][i][j][k] : flux
	 * localNE[2][i][j][k] : conc
	 * localNE[3][i][j][k] : frac
	 */
	// virtual std::vector<std::vector<std::vector<std::tuple<double, double, double, double> > > > & getLocalNE(
	// 		std::shared_ptr<xolotlSolver::PetscSolver> solver) {
	// 	//std::tuple<double, double, double, double> v = {0.0,0.0,0.0,0.0};
	// 	std::vector<std::vector<std::vector<std::tuple<double, double, double, double> > > > a;
	// 	return a;
	// };

	// /**
	//  * Get the local Xe rate that needs to be passed
	//  *
	//  * @param solver The pointer to the solver
	//  * @return The local vector of rates
	//  */
	// virtual std::vector<std::vector<std::vector<double> > > * getLocalXeRate(
	// 		std::shared_ptr<xolotlSolver::PetscSolver> solver) {
	// 	std::vector<std::vector<std::vector<double> > > *a = {NULL};
	// 	return a;
	// };


	virtual double getLocalXeRatePoint(
			std::shared_ptr<xolotlSolver::PetscSolver> solver,
			int i, int j, int k) {
		return 0.0;
	};

	virtual double getLocalXeFluxPoint(
			std::shared_ptr<xolotlSolver::PetscSolver> solver,
			int i, int j, int k) {
		return 0.0;
	};

	virtual double getLocalXeConcPoint(
			std::shared_ptr<xolotlSolver::PetscSolver> solver,
			int i, int j, int k) {
		return 0.0;
	};

	virtual double getLocalXeVolFracPoint(
			std::shared_ptr<xolotlSolver::PetscSolver> solver,
			int i, int j, int k) {
		return 0.0;
	};


	/**
	 * Get the local concentration that needs to be passed
	 *
	 * @param solver The pointer to the solver
	 * @return The local vector of concentrations
	 */
	// virtual std::vector<std::vector<std::vector<double> > > * getLocalXeConc(
	// 		std::shared_ptr<xolotlSolver::PetscSolver> solver) {
	// 	std::vector<std::vector<std::vector<double> > > *a = {NULL};
	// 	return a;
	// };

	/**
	 * Get the local coordinate that needs to be passed
	 *
	 * @param solver The pointer to the solver
	 * @param xs, xm The start and width in the X direction on the local MPI process
	 * @param Mx The total width in the X direction
	 * @param ys, ym The start and width in the Y direction on the local MPI process
	 * @param My The total width in the Y direction
	 * @param zs, zm The start and width in the Z direction on the local MPI process
	 * @param Mz The total width in the Z direction
	 */
	virtual void getLocalCoordinates(std::shared_ptr<xolotlSolver::PetscSolver> solver,
	 		int *xs, int *xm, int *Mx, int *ys, int *ym, int *My, int *zs,
	 		int *zm, int *Mz){};

	virtual void initGBLocation(
		std::shared_ptr<xolotlSolver::PetscSolver> solver){};

	// virtual void setGBLocations(
	// 		std::shared_ptr<xolotlSolver::PetscSolver> solver,
	// 		std::vector<std::tuple<int, int, int> > gbList) {};

	virtual void setGBLocation(
			std::shared_ptr<xolotlSolver::PetscSolver> solver, int i, int j,
			int k) {};

	/**
	 * Reset the GB vector.
	 */
	virtual void resetGBVector(std::shared_ptr<xolotlSolver::PetscSolver> solver){};

	/**
	 * Finalize the solve
	 */
	virtual void finalizeXolotl(std::shared_ptr<xolotlSolver::PetscSolver> solver, bool isStandalone = true){};

 /**
	* Get the TS from the solver
	*
	* @return The TS
	*/
	// virtual TS getXolotlTS(std::shared_ptr<xolotlSolver::PetscSolver> solver) {
  // 	TS toReturn = nullptr;
	// 	return toReturn;
  // };

};
// End class interface

#endif
