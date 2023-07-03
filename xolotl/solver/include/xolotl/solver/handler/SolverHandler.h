#pragma once

// Includes
#include <xolotl/core/Constants.h>
#include <xolotl/solver/handler/ISolverHandler.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/RandomNumberGenerator.h>

namespace xolotl
{
namespace solver
{
namespace handler
{
/**
 * This class and its subclasses realize the ISolverHandler interface to solve
 * the advection-diffusion-reaction problem with currently supported solvers.
 */
class SolverHandler : public ISolverHandler
{
public:
	using NetworkType = core::network::IReactionNetwork;

	using ConcentrationsView = typename NetworkType::ConcentrationsView;
	using FluxesView = typename NetworkType::FluxesView;
	using SparseFillMap = typename NetworkType::SparseFillMap;

protected:
	/**
	 * The vector to know where the GB are.
	 *
	 * The first pair is the location of a grid point (X,Y),
	 */
	std::vector<std::array<IdType, 3>> gbVector;

	//! The name of the network file
	std::string networkName;

	//! The name of the free GB file
	std::string gbFileName;

	//! The original network created from the network loader.
	NetworkType& network;

	//! Vector storing the grid in the x direction
	std::vector<double> grid;

	//! Vector storing the previous grid in the x direction
	std::vector<double> oldGrid;

	//! Vector storing the grid for the temperature
	std::vector<double> temperatureGrid;

	//! The current temperature on the grid.
	std::vector<double> temperature;

	//! The number of grid points in the depth direction.
	IdType nX;

	//! The number of grid points in the Y direction.
	IdType nY;

	//! The number of grid points in the Z direction.
	IdType nZ;

	//! The grid step size in the y direction.
	double hY;

	//! The grid step size in the z direction.
	double hZ;

	//! The power used for the temperature grid step size.
	double tempGridPower;

	//! The local start of grid points in the X direction.
	IdType localXS;

	//! The local width of grid points in the X direction.
	IdType localXM;

	//! The local start of grid points in the Y direction.
	IdType localYS;

	//! The local width of grid points in the Y direction.
	IdType localYM;

	//! The local start of grid points in the Z direction.
	IdType localZS;

	//! The local width of grid points in the Z direction.
	IdType localZM;

	//! The number of grid points by which the boundary condition should be
	//! shifted at this side.
	IdType leftOffset, rightOffset, bottomOffset, topOffset, frontOffset,
		backOffset;

	//! The initial vacancy concentration.
	std::vector<std::pair<IdType, double>> initialConc;

	//! The vector of quantities to pass to MOOSE.
	// 0: Xe rate, 1: previous flux, 2: monomer concentration, 3: volume
	// fraction
	std::vector<std::vector<std::vector<std::array<double, 4>>>> localNE;

	//! The electronic stopping power for re-solution
	double electronicStoppingPower;

	//! The original flux handler created.
	core::flux::IFluxHandler* fluxHandler;

	//! The original temperature handler created.
	core::temperature::ITemperatureHandler* temperatureHandler;

	//! the original perf handler created.
	std::shared_ptr<perf::IPerfHandler> perfHandler;

	//! the original viz handler created.
	std::shared_ptr<viz::IVizHandler> vizHandler;

	//! The original diffusion handler created.
	core::diffusion::IDiffusionHandler* diffusionHandler;

	//! The original Soret diffusion handler created.
	core::modified::ISoretDiffusionHandler* soretDiffusionHandler;

	//! The vector of advection handlers.
	std::vector<core::advection::IAdvectionHandler*> advectionHandlers;

	//! The number of dimensions for the problem.
	int dimension;

	//! If the user wants to move the surface.
	bool movingSurface;

	//! If the user wants to use x mirror boundary conditions or periodic ones.
	bool isMirror;

	//! If the user wants to use x Robin boundary conditions for temperature.
	bool isRobin;

	//! If the user wants to attenuate the modified trap mutation.
	bool useAttenuation;

	//! What type of temperature grid to use.
	bool sameTemperatureGrid;

	//! If the user wants to use a temporal profile for the flux.
	bool fluxTempProfile;

	//! The sputtering yield for the problem.
	double sputteringYield;

	//! The depth parameter for the bubble bursting.
	double tauBursting;

	//! The factor involved in computing bursting likelihood.
	double burstingFactor;

	//! The value to use to seed the random number generator.
	unsigned int rngSeed;

	//! The minimum sizes for average radius computation.
	std::vector<size_t> minRadiusSizes;

	//! The previous time.
	double previousTime;

	//! The number of xenon atoms that went to the GB
	double nXeGB;

	//! The grid options
	std::string gridType;
	std::string gridFileName;
	double gridParam0, gridParam1, gridParam2, gridParam3, gridParam4,
		gridParam5;

	//! The random number generator to use.
	std::unique_ptr<util::RandomNumberGenerator<int, unsigned int>> rng;

	/**
	 * Method generating the grid in the x direction
	 *
	 * @param surfaceOffset The number of grid point to add/remove at the
	 * surface
	 */
	void
	generateGrid(int surfaceOffset);

	/**
	 * Constructor.
	 *
	 * @param _network The reaction network to use.
	 * @param options The options.
	 */
	SolverHandler(NetworkType& _network, const options::IOptions& options);

public:
	//! The Constructor
	SolverHandler() = delete;

	virtual ~SolverHandler();

	/**
	 * \see ISolverHandler.h
	 */
	void
	initializeHandlers(core::material::IMaterialHandler* material,
		core::temperature::ITemperatureHandler* tempHandler,
		const options::IOptions& opts) override;

	/**
	 * \see ISolverHandler.h
	 */
	void
	generateTemperatureGrid() override
	{
		// Don't do anything if we want the same grid as the cluster one
		if (sameTemperatureGrid) {
			temperatureGrid = grid;
			return;
		}

		// If the temperature grid already existed we need to save its values
		std::vector<double> oldGrid;
		if (temperatureGrid.size() > 0)
			oldGrid = temperatureGrid;
		// Clear the grid
		temperatureGrid.clear();

		// Doesn't mean anything in 0D
		if (grid.size() == 0)
			return;

		// Compute the total width
		auto n = grid.size() - 2;
		auto width =
			((grid[grid.size() - 3] + grid[grid.size() - 2]) / 2.0 - grid[1]);

		auto newWidth = width;
		auto newH = pow(newWidth, 1 / tempGridPower) / (n - 1.5);

		// Surface
		temperatureGrid.push_back(0);
		temperatureGrid.push_back(pow(newH, tempGridPower));

		// Material
		for (auto i = 2; i < grid.size(); i++) {
			auto j = i - 1;
			temperatureGrid.push_back(
				temperatureGrid[1] + (pow(j * newH, tempGridPower)));
		}

		// The temperature values need to be updated to match the new grid
		if (oldGrid.size() == 0)
			return;

		// Get the default temperature vector if needed
		auto localTemp = temperature;
		// First, broadcast the temperature vector so that each
		// rank can access any temperature on the grid.
		auto xolotlComm = util::getMPIComm();
		int procId;
		MPI_Comm_rank(xolotlComm, &procId);
		int worldSize;
		MPI_Comm_size(xolotlComm, &worldSize);

		// Need to create an array of size worldSize where each entry
		// is the number of element on the given procId
		int toSend = localXM;
		// Ghost cells
		if (localXS == 0 || localXS + localXM == nX)
			toSend++;
		if (localXS == 0 && localXS + localXM == nX)
			toSend++;

		// Receiving array for number of elements
		int counts[worldSize];
		MPI_Allgather(&toSend, 1, MPI_INT, counts, 1, MPI_INT, xolotlComm);

		// Define the displacements
		int displacements[worldSize];
		for (auto i = 0; i < worldSize; i++) {
			if (i == 0)
				displacements[i] = 0;
			else
				displacements[i] = displacements[i - 1] + counts[i - 1];
		}

		// Receiving array for temperatures
		double broadcastedTemp[nX + 2];
		double valuesToSend[toSend];
		if (localXS == 0) {
			for (auto i = 0; i < toSend; i++) {
				valuesToSend[i] = localTemp[i];
			}
		}
		else {
			for (auto i = 0; i < toSend; i++) {
				valuesToSend[i] = localTemp[i + 1];
			}
		}
		MPI_Allgatherv(&valuesToSend, toSend, MPI_DOUBLE, &broadcastedTemp,
			counts, displacements, MPI_DOUBLE, xolotlComm);

		// Now we can interpolate the temperature for the local grid
		std::vector<double> toReturn;
		// Loop on the local grid including ghosts
		for (auto i = localXS; i < localXS + localXM + 2; i++) {
			// Left of surface
			if (i < 0) {
				toReturn.push_back(broadcastedTemp[nX]);
				continue;
			}
			// Get the grid location
			double loc = 0.0;
			if (i == 0)
				loc = temperatureGrid[0] - temperatureGrid[1];
			else
				loc = (temperatureGrid[i - 1] + temperatureGrid[i]) / 2.0 -
					temperatureGrid[1];

			bool matched = false;
			IdType jKeep = 0;
			// Look for it in the temperature grid
			for (auto j = jKeep; j < nX + 1; j++) {
				double tempLoc1 = 0.0,
					   tempLoc2 =
						   (oldGrid[j] + oldGrid[j + 1]) / 2.0 - oldGrid[1];
				if (j == 0)
					tempLoc1 = oldGrid[0] - oldGrid[1];
				else
					tempLoc1 = (oldGrid[j - 1] + oldGrid[j]) / 2.0 - oldGrid[1];

				if (loc >= tempLoc1 && loc < tempLoc2) {
					double xLoc = (loc - tempLoc1) / (tempLoc2 - tempLoc1);
					double y1 = broadcastedTemp[j], y2 = broadcastedTemp[j + 1];
					toReturn.push_back(y1 + xLoc * (y2 - y1));
					matched = true;
					jKeep = j;
					break;
				}
			}

			if (not matched)
				toReturn.push_back(broadcastedTemp[nX]);
		}

		temperature = toReturn;
	}

	/**
	 * \see ISolverHandler.h
	 */
	std::vector<double>
	getXGrid() const override
	{
		return grid;
	}

	/**
	 * \see ISolverHandler.h
	 */
	std::vector<double>
	getTemperatureGrid() const override
	{
		return temperatureGrid;
	}

	/**
	 * \see ISolverHandler.h
	 */
	double
	getStepSizeY() const override
	{
		return hY;
	}

	/**
	 * \see ISolverHandler.h
	 */
	double
	getStepSizeZ() const override
	{
		return hZ;
	}

	/**
	 * \see ISolverHandler.h
	 */
	int
	getDimension() const override
	{
		return dimension;
	}

	/**
	 * \see ISolverHandler.h
	 */
	std::vector<std::pair<IdType, double>>
	getInitialConc() const override
	{
		return initialConc;
	}

	/**
	 * \see ISolverHandler.h
	 */
	double
	getSputteringYield() const override
	{
		return sputteringYield;
	}

	/**
	 * \see ISolverHandler.h
	 */
	double
	getTauBursting() const override
	{
		return tauBursting;
	}

	/**
	 * \see ISolverHandler.h
	 */
	double
	getBurstingFactor() const override
	{
		return burstingFactor;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	createLocalNE(IdType a, IdType b = 1, IdType c = 1) override;

	/**
	 * \see ISolverHandler.h
	 */
	void
	setLocalXeRate(double rate, IdType i, IdType j = 0, IdType k = 0) override
	{
		std::get<0>(localNE[i][j][k]) += rate;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setLocalNE(
		const std::vector<std::vector<std::vector<std::array<double, 4>>>>&
			rateVector) override
	{
		localNE = rateVector;
	}

	/**
	 * \see ISolverHandler.h
	 */
	std::vector<std::vector<std::vector<std::array<double, 4>>>>&
	getLocalNE() override
	{
		return localNE;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setPreviousXeFlux(
		double flux, IdType i, IdType j = 0, IdType k = 0) override
	{
		std::get<1>(localNE[i][j][k]) = flux;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setMonomerConc(double conc, IdType i, IdType j = 0, IdType k = 0) override
	{
		std::get<2>(localNE[i][j][k]) = conc;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setVolumeFraction(
		double frac, IdType i, IdType j = 0, IdType k = 0) override
	{
		std::get<3>(localNE[i][j][k]) = frac;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setLocalCoordinates(IdType xs, IdType xm, IdType ys = 0, IdType ym = 0,
		IdType zs = 0, IdType zm = 0) override
	{
		localXS = xs;
		localXM = xm;
		localYS = ys;
		localYM = ym;
		localZS = zs;
		localZM = zm;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	getLocalCoordinates(IdType& xs, IdType& xm, IdType& Mx, IdType& ys,
		IdType& ym, IdType& My, IdType& zs, IdType& zm, IdType& Mz) override
	{
		xs = localXS;
		xm = localXM;
		Mx = nX;
		ys = localYS;
		ym = localYM;
		My = nY;
		zs = localZS;
		zm = localZM;
		Mz = nZ;
	}

	/**
	 * \see ISolverHandler.h
	 */
	IdType
	getLeftOffset() const override
	{
		return leftOffset;
	}

	/**
	 * \see ISolverHandler.h
	 */
	IdType
	getRightOffset() const override
	{
		return rightOffset;
	}

	/**
	 * \see ISolverHandler.h
	 */
	bool
	moveSurface() const override
	{
		return movingSurface;
	}

	/**
	 * \see ISolverHandler.h
	 */
	bool
	temporalFlux() const override
	{
		return fluxTempProfile;
	}

	/**
	 * \see ISolverHandler.h
	 */
	double
	getPreviousTime() override
	{
		return previousTime;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setPreviousTime(double time, bool updateFluence = false) override
	{
		previousTime = time;
		if (updateFluence)
			fluxHandler->computeFluence(time);
	}

	/**
	 * \see ISolverHandler.h
	 */
	double
	getNXeGB() override
	{
		return nXeGB;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setNXeGB(double nXe) override
	{
		nXeGB = nXe;
	}

	/**
	 * \see ISolverHandler.h
	 */
	std::vector<size_t>
	getMinSizes() const override
	{
		return minRadiusSizes;
	}

	/**
	 * \see ISolverHandler.h
	 */
	core::flux::IFluxHandler*
	getFluxHandler() const override
	{
		return fluxHandler;
	}

	/**
	 * \see ISolverHandler.h
	 */
	core::temperature::ITemperatureHandler*
	getTemperatureHandler() const override
	{
		return temperatureHandler;
	}

	/**
	 * \see ISolverHandler.h
	 */
	std::shared_ptr<perf::IPerfHandler>
	getPerfHandler() const override
	{
		return perfHandler;
	}

	/**
	 * \see ISolverHandler.h
	 */
	std::shared_ptr<viz::IVizHandler>
	getVizHandler() const override
	{
		return vizHandler;
	}

	/**
	 * \see ISolverHandler.h
	 */
	core::diffusion::IDiffusionHandler*
	getDiffusionHandler() const override
	{
		return diffusionHandler;
	}

	/**
	 * \see ISolverHandler.h
	 */
	core::modified::ISoretDiffusionHandler*
	getSoretDiffusionHandler() const override
	{
		return soretDiffusionHandler;
	}

	/**
	 * \see ISolverHandler.h
	 */
	core::advection::IAdvectionHandler*
	getAdvectionHandler() const override
	{
		return advectionHandlers[0];
	}

	/**
	 * \see ISolverHandler.h
	 */
	std::vector<core::advection::IAdvectionHandler*>
	getAdvectionHandlers() const override
	{
		return advectionHandlers;
	}

	/**
	 * \see ISolverHandler.h
	 */
	virtual core::network::IReactionNetwork&
	getNetwork() const override
	{
		return network;
	}

	/**
	 * \see ISolverHandler.h
	 */
	std::string
	getNetworkName() const override
	{
		return networkName;
	}

	/**
	 * \see ISolverHandler.h
	 */
	util::RandomNumberGenerator<int, unsigned int>&
	getRNG(void) const override
	{
		return *rng;
	}

	/**
	 * \see ISolverHandler.h
	 */
	std::vector<std::array<IdType, 3>>
	getGBVector() const override
	{
		return gbVector;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setGBLocation(IdType i, IdType j = 0, IdType k = 0) override;

	/**
	 * \see ISolverHandler.h
	 */
	void
	resetGBVector() override
	{
		gbVector.clear();
	}

	/**
	 * \see ISolverHandler.h
	 */
	std::vector<double>
	interpolateTemperature(
		std::vector<double> localTemp = std::vector<double>()) override
	{
		// No need to interpolate if the grid are the same
		if (sameTemperatureGrid)
			return temperature;

		// Get the default temperature vector if needed
		if (localTemp.size() == 0)
			localTemp = temperature;
		// First, broadcast the temperature vector so that each
		// rank can access any temperature on the grid.
		auto xolotlComm = util::getMPIComm();
		int procId;
		MPI_Comm_rank(xolotlComm, &procId);
		int worldSize;
		MPI_Comm_size(xolotlComm, &worldSize);

		// Need to create an array of size worldSize where each entry
		// is the number of element on the given procId
		int toSend = localXM;
		// Ghost cells
		if (localXS == 0 || localXS + localXM == nX)
			toSend++;
		if (localXS == 0 && localXS + localXM == nX)
			toSend++;

		// Receiving array for number of elements
		int counts[worldSize];
		MPI_Allgather(&toSend, 1, MPI_INT, counts, 1, MPI_INT, xolotlComm);

		// Define the displacements
		int displacements[worldSize];
		for (auto i = 0; i < worldSize; i++) {
			if (i == 0)
				displacements[i] = 0;
			else
				displacements[i] = displacements[i - 1] + counts[i - 1];
		}

		// Receiving array for temperatures
		double broadcastedTemp[nX + 2];
		double valuesToSend[toSend];
		if (localXS == 0) {
			for (auto i = 0; i < toSend; i++) {
				valuesToSend[i] = localTemp[i];
			}
		}
		else {
			for (auto i = 0; i < toSend; i++) {
				valuesToSend[i] = localTemp[i + 1];
			}
		}
		MPI_Allgatherv(&valuesToSend, toSend, MPI_DOUBLE, &broadcastedTemp,
			counts, displacements, MPI_DOUBLE, xolotlComm);

		// Now we can interpolate the temperature for the local grid
		std::vector<double> toReturn;
		// Loop on the local grid including ghosts
		for (auto i = localXS; i < localXS + localXM + 2; i++) {
			// Get the grid location
			double loc = 0.0;
			if (i == 0)
				loc = grid[0] - grid[1];
			else
				loc = (grid[i - 1] + grid[i]) / 2.0 - grid[1];

			bool matched = false;
			IdType jKeep = 0;
			// Look for it in the temperature grid
			for (auto j = jKeep; j < nX + 1; j++) {
				double tempLoc1 = 0.0,
					   tempLoc2 =
						   (temperatureGrid[j] + temperatureGrid[j + 1]) / 2.0 -
					temperatureGrid[1];
				if (j == 0)
					tempLoc1 = temperatureGrid[0] - temperatureGrid[1];
				else
					tempLoc1 =
						(temperatureGrid[j - 1] + temperatureGrid[j]) / 2.0 -
						temperatureGrid[1];

				if (loc >= tempLoc1 && loc < tempLoc2) {
					double xLoc = (loc - tempLoc1) / (tempLoc2 - tempLoc1);
					double y1 = broadcastedTemp[j], y2 = broadcastedTemp[j + 1];
					toReturn.push_back(y1 + xLoc * (y2 - y1));
					matched = true;
					jKeep = j;
					break;
				}
			}

			if (not matched) {
				toReturn.push_back(broadcastedTemp[nX + 1]);
			}
		}

		return toReturn;
	}
};
} /* namespace handler */
} /* namespace solver */
} /* namespace xolotl */
