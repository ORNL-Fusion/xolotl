#pragma once

#include <petscts.h>

#include <memory>
#include <vector>

#include <mpi.h>

#include <xolotl/config.h>
#include <xolotl/interface/IXolotlInterface.h>

namespace xolotl
{
namespace options
{
class IOptions;
}
namespace core
{
namespace network
{
class IReactionNetwork;
}
} // namespace core
namespace solver
{
class ISolver;
}
namespace perf
{
class IPerfHandler;
}

namespace interface
{
class ComputeContext;

/**
 * Class defining the method to be coupled to another code through MOOSEApps
 */
class XolotlInterface : public IXolotlInterface
{
private:
	/**
	 * Did this object initialize xolotl?
	 */
	bool initializedHere{false};

	/**
	 * Did this instance initialize the solver?
	 */
	bool solverInitialized{false};

	/**
	 * The MPI and Kokkos environment
	 */
	std::shared_ptr<ComputeContext> computeContext;

	/**
	 * The solver
	 */
	std::shared_ptr<solver::ISolver> solver;

	/**
	 * A vector of maps to know which cluster in subnetworks correspond
	 * to which in the main network
	 */
	std::vector<std::vector<IdType>> fromSubNetwork;

	/**
	 * A (possibly) shared instance of the options object
	 */
	std::shared_ptr<options::IOptions> options;

public:
	/**
	 * The default constructor
	 */
	XolotlInterface();

	/**
	 * Construct from command-line arguments
	 *
	 * @param argc, argv The command line arguments
	 * @param mpiComm The communicator to use
	 */
	XolotlInterface(
		int& argc, const char* argv[], MPI_Comm mpiComm = MPI_COMM_WORLD);

	/**
	 * Construct from preconfigured compute context and options
	 *
	 * @param opts The options object
	 * @param mpiComm The communicator to use
	 */
	XolotlInterface(const std::shared_ptr<ComputeContext>& context,
		const std::shared_ptr<options::IOptions>& opts,
		MPI_Comm mpiComm = MPI_COMM_WORLD);

	/**
	 * The destructor
	 */
	virtual ~XolotlInterface();

	/**
	 * Print something
	 */
	void
	printSomething();

	/**
	 * Initialize all the options and handlers
	 *
	 * @param argc, argv The command line arguments
	 * @param comm The communicator to use
	 */
	void
	initializeXolotl(
		int& argc, const char* argv[], MPI_Comm comm = MPI_COMM_WORLD);

	/**
	 * Initialize all the options and handlers
	 *
	 * @note Assumes options and MPI communicator have already been set
	 */
	void
	initializeXolotl();

	void
	initializeSolver();

	/**
	 * Set the final time and the dt.
	 *
	 * @param finalTime The wanted final time
	 * @param dt The wanted max time step
	 */
	void
	setTimes(double finalTime, double dt);

	/**
	 * Run the PETSc solve
	 */
	void
	solveXolotl() override;

	/**
	 * Get the vector of data that can be passed to an app
	 *
	 * @return The vector
	 */
	std::vector<std::vector<std::vector<std::array<double, 4>>>>
	getLocalNE();

	/**
	 * Set the vector of data from an app
	 *
	 * @param The vector
	 */
	void
	setLocalNE(
		const std::vector<std::vector<std::vector<std::array<double, 4>>>>&
			rateVector);

	/**
	 * Get the local coordinates of the grid, sizes and indices
	 *
	 * @param xs, xm The start and width in the X direction on the local MPI
	 * process
	 * @param Mx The total width in the X direction
	 * @param ys, ym The start and width in the Y direction on the local MPI
	 * process
	 * @param My The total width in the Y direction
	 * @param zs, zm The start and width in the Z direction on the local MPI
	 * process
	 * @param Mz The total width in the Z direction
	 */
	void
	getLocalCoordinates(IdType& xs, IdType& xm, IdType& Mx, IdType& ys,
		IdType& ym, IdType& My, IdType& zs, IdType& zm, IdType& Mz);

	/**
	 * Set the location of one GB grid point.
	 *
	 * @param i, j, k The coordinate of the GB
	 */
	void
	setGBLocation(IdType i, IdType j = 0, IdType k = 0);

	/**
	 * Reset the GB vector.
	 */
	void
	resetGBVector();

	/**
	 * Get the concentrations and their ids.
	 *
	 * @return The concentration vector from the current state of the simulation
	 */
	std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
	getConcVector();

	/**
	 * Set the concentrations and their ids.
	 *
	 * @ param concVector A given state of the concentrations
	 */
	void
	setConcVector(std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
			concVector);

	/**
	 * Get the previous time.
	 *
	 * @return The previous time
	 */
	double
	getPreviousTime();

	/**
	 * Set the previous time.
	 *
	 * @param time The previous time
	 */
	void
	setPreviousTime(double time);

	/**
	 * Get the current dt.
	 *
	 * @return  The current time step
	 */
	double
	getCurrentDt();

	/**
	 * Set the current time and dt.
	 *
	 * @param currentTime The time
	 * @param currentDt The current time step
	 */
	void
	setCurrentTimes(double currentTime, double currentDt);

	/**
	 * Get the number of Xe that went to the GB.
	 *
	 * @return The number of Xenon
	 */
	double
	getNXeGB();

	/**
	 * Set the number of Xe that went to the GB.
	 *
	 * @param nXe The number of Xenon
	 */
	void
	setNXeGB(double nXe);

	/**
	 * Get the TS from the solver.
	 *
	 * @return The TS
	 */
	TS&
	getTS();

	/**
	 * Get the grid information
	 *
	 * @param hy The spacing in the Y direction
	 * @param hz The spacing in the Z direction
	 * @return The grid in the X direction
	 */
	std::vector<double>
	getGridInfo(double& hy, double& hz);

	/**
	 * Get the cluster information
	 *
	 * @return The vector representing the bounds for each cluster
	 */
	std::vector<std::vector<AmountType>>
	getAllClusterBounds();

	/**
	 * Get the cluster information
	 *
	 * @return The vector linking moment Ids to cluster Ids
	 */
	std::vector<std::vector<IdType>>
	getAllMomentIdInfo();

	/**
	 * Computes the map between the different set of cluster bounds and moment
	 * IDs.
	 *
	 * @param bounds A vector of cluster bounds
	 * @param momIdInfo Information about which moment ID goes with which
	 * cluster
	 */
	void
	initializeClusterMaps(
		std::vector<std::vector<std::vector<AmountType>>> bounds,
		std::vector<std::vector<std::vector<IdType>>> momIdInfo);

	/**
	 * Initializes the reaction in a separate step.
	 */
	void
	initializeReactions();

	/**
	 * Get the implanted flux for each sub network.
	 *
	 * @return The vector of vectors of flux, first is the ID and second is the
	 * value.
	 */
	std::vector<std::vector<std::pair<IdType, double>>>
	getImplantedFlux();

	/**
	 * Set the implanted flux for each sub network.
	 *
	 * @param fluxVector With first is the ID and second is the value
	 */
	void
	setImplantedFlux(std::vector<std::pair<IdType, double>> fluxVector);

	/**
	 * Values for the rates to be set in constant reactions.
	 *
	 * @param rates All the rates
	 */
	void
	setConstantRates(std::vector<std::vector<double>> rates);

	/**
	 * Compute the constant rates
	 *
	 * @param conc The concentration vector
	 * @return A vector containing the rates for each sub instance
	 */
	std::vector<std::vector<std::vector<double>>>
	computeConstantRates(std::vector<std::vector<double>> conc);

	/**
	 * Get the connectivity matrices
	 *
	 * @return A vector telling which reactants interact together
	 */
	std::vector<std::vector<std::vector<bool>>>
	getConstantConnectivities();

	/**
	 * Set the connectivity matrices for constant reactions
	 *
	 * @param conns A vector telling which reactants interact together
	 */
	void
	setConstantConnectivities(std::vector<std::vector<bool>> conns);

	/**
	 * Write the data in a file.
	 *
	 * @param time The current time
	 * @param conc The concentration vector
	 */
	void
	outputData(double time, std::vector<std::vector<double>> conc);

	/**
	 * Get whether the solve converged or not.
	 *
	 * @return true if it converged
	 */
	bool
	getConvergenceStatus();

	/**
	 * Finalize the solve
	 */
	void
	finalizeXolotl();
};
} /* namespace interface */
} /* namespace xolotl */
