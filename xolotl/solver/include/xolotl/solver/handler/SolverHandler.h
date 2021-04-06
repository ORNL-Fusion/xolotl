#ifndef SOLVERHANDLER_H
#define SOLVERHANDLER_H

// Includes
#include <xolotl/core/Constants.h>
#include <xolotl/io/XFile.h>
#include <xolotl/solver/handler/ISolverHandler.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/RandomNumberGenerator.h>
#include <xolotl/util/TokenizedLineReader.h>

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
	double initialVConc;

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

	//! The original diffusion handler created.
	core::diffusion::IDiffusionHandler* diffusionHandler;

	//! The vector of advection handlers.
	std::vector<core::advection::IAdvectionHandler*> advectionHandlers;

	//! The original modified trap-mutation handler created.
	core::modified::ITrapMutationHandler* mutationHandler;

	//! The number of dimensions for the problem.
	int dimension;

	//! The portion of void at the beginning of the problem.
	double portion;

	//! If the user wants to move the surface.
	bool movingSurface;

	//! If the user wants to burst bubbles.
	bool bubbleBursting;

	//! If the user wants to use x mirror boundary conditions or periodic ones.
	bool isMirror;

	//! If the user wants to attenuate the modified trap mutation.
	bool useAttenuation;

	//! The sputtering yield for the problem.
	double sputteringYield;

	//! The depth parameter for the bubble bursting.
	double tauBursting;

	//! The factor involved in computing bursting likelihood.
	double burstingFactor;

	//! The ratio of He per V in a bubble.
	double heVRatio;

	//! The value to use to seed the random number generator.
	unsigned int rngSeed;

	//! The minimum sizes for average radius computation.
	std::vector<size_t> minRadiusSizes;

	//! The previous time.
	double previousTime;

	//! The number of xenon atoms that went to the GB
	double nXeGB;

	//! The random number generator to use.
	std::unique_ptr<util::RandomNumberGenerator<int, unsigned int>> rng;

	/**
	 * Method generating the grid in the x direction
	 *
	 * @param opts The options
	 */
	void
	generateGrid(const options::IOptions& opts)
	{
		// Clear the grid
		grid.clear();

		// Doesn't mean anything in 0D
		if (opts.getDimensionNumber() == 0)
			return;

		// Get the grid type
		auto gridType = opts.getGridTypeName();

		// Check if we want to read in the grid from a file
		if (gridType == "read") {
			// Open the corresponding file
			std::ifstream inputFile(opts.getGridFilename().c_str());
			if (!inputFile)
				std::cerr << "\nCould not open the file containing the grid "
							 "spacing information. "
							 "Aborting!\n"
						  << std::endl;
			// Get the data
			std::string line;
			getline(inputFile, line);

			// Break the line into a vector
			util::TokenizedLineReader<double> reader;
			auto argSS = std::make_shared<std::istringstream>(line);
			reader.setInputStream(argSS);
			auto tokens = reader.loadLine();

			if (tokens.size() == 0)
				std::cerr << "\nDid not read correctly the file containing the "
							 "grid spacing information. "
							 "Aborting!\n"
						  << std::endl;

			// Compute the offset to add to the grid for boundary conditions
			double offset = tokens[1] - tokens[0];
			// Add the first grid point
			grid.push_back(0.0);
			// Check the location of the first grid point
			if (tokens[0] > 0.0) {
				grid.push_back(offset);
			}

			// Loop on the tokens
			for (auto i = 0; i < tokens.size(); i++) {
				grid.push_back(tokens[i] + offset);
			}

			// Add the last grid point for boundary conditions
			grid.push_back(2.0 * tokens[tokens.size() - 1] -
				tokens[tokens.size() - 2] + offset);

			// Set the number of grid points
			nX = grid.size() - 2;

			// Get the number of dimensions
			auto dim = opts.getDimensionNumber();
			if (dim > 1) {
				nY = opts.getGridParam(0);
				hY = opts.getGridParam(1);
			}
			if (dim > 2) {
				nZ = opts.getGridParam(2);
				hZ = opts.getGridParam(3);
			}

			return;
		}

		// Maybe the user wants a Chebyshev grid
		if (gridType == "cheby") {
			// The first grid point will be at x = 0.0
			grid.push_back(0.0);
			grid.push_back(0.0);

			// Set the number of grid points
			nX = opts.getGridParam(0);
			auto hx = opts.getGridParam(1);
			// In that case hx correspond to the full length of the grid
			for (auto l = 1; l <= nX - 1; l++) {
				grid.push_back((hx / 2.0) *
					(1.0 - cos(core::pi * double(l) / double(nX - 1))));
			}
			// The last grid point will be at x = hx
			grid.push_back(hx);

			// Get the number of dimensions
			auto dim = opts.getDimensionNumber();
			if (dim > 1) {
				nY = opts.getGridParam(2);
				hY = opts.getGridParam(3);
			}
			if (dim > 2) {
				nZ = opts.getGridParam(4);
				hZ = opts.getGridParam(5);
			}

			return;
		}
		// Check if the user wants a regular grid
		if (gridType == "uniform") {
			// Set the number of grid points
			nX = opts.getGridParam(0);
			auto hx = opts.getGridParam(1);
			// The grid will me made of nx + 2 points separated by hx nm
			for (auto l = 0; l <= nX + 1; l++) {
				grid.push_back((double)l * hx);
			}

			// Get the number of dimensions
			auto dim = opts.getDimensionNumber();
			if (dim > 1) {
				nY = opts.getGridParam(2);
				hY = opts.getGridParam(3);
			}
			if (dim > 2) {
				nZ = opts.getGridParam(4);
				hZ = opts.getGridParam(5);
			}

			return;
		}

		// If it is not regular do a fine mesh close to the surface and
		// increase the step size when away from the surface
		if (gridType == "nonuniform") {
			// Initialize the value of the previous point
			double previousPoint = 0.0;
			// Set the number of grid points
			nX = opts.getGridParam(0);
			// Set the position of the surface
			auto surfacePos = (IdType)(nX * opts.getVoidPortion() / 100.0);

			// Loop on all the grid points
			for (auto l = 0; l <= nX + 1; l++) {
				// Add the previous point
				grid.push_back(previousPoint);
				// 0.1nm step near the surface (x < 2.5nm)
				if (l < surfacePos + 26) {
					previousPoint += 0.1;
				}
				// Then 0.25nm (2.5nm < x < 5.0nm)
				else if (l < surfacePos + 36) {
					previousPoint += 0.25;
				}
				// Then 0.5nm (5.0nm < x < 7.5nm)
				else if (l < surfacePos + 41) {
					previousPoint += 0.5;
				}
				// Then 1.0nm step size (7.5nm < x < 50.5)
				else if (l < surfacePos + 84) {
					previousPoint += 1.0;
				}
				// Then 2.0nm step size (50.5nm < x < 100.5)
				else if (l < surfacePos + 109) {
					previousPoint += 2.0;
				}
				// Then 5.0nm step size (100.5nm < x < 150.5)
				else if (l < surfacePos + 119) {
					previousPoint += 5.0;
				}
				// Then 10.0nm step size (150.5nm < x < 300.5)
				else if (l < surfacePos + 134) {
					previousPoint += 10.0;
				}
				// Then 20.0nm step size (300.5nm < x < 500.5)
				else if (l < surfacePos + 144) {
					previousPoint += 20.0;
				}
				// Then 50.0nm step size (500.5nm < x < 1000.5)
				else if (l < surfacePos + 154) {
					previousPoint += 50.0;
				}
				// Then 100.0nm step size (1000.5nm < x < 5000.5)
				else if (l < surfacePos + 194) {
					previousPoint += 100.0;
				}
				// Then 200.0nm step size (5000.5nm < x < 10000.5)
				else if (l < surfacePos + 219) {
					previousPoint += 200.0;
				}
				// Then 500.0nm step size (10000.5nm < x < 20000.5)
				else if (l < surfacePos + 239) {
					previousPoint += 500.0;
				}
				// Then 1.0um step size (20000.5nm < x < 30000.5nm )
				else if (l < surfacePos + 249) {
					previousPoint += 1000.0;
				}
				// Then 2.0um step size (30000.5nm < x < 50000.5)
				else if (l < surfacePos + 259) {
					previousPoint += 2000.0;
				}
				// Then 5.0um step size (50000.5nm < x < 100000.5)
				else if (l < surfacePos + 269) {
					previousPoint += 5000.0;
				}
				// Then 10.0um step size (100000.5nm < x < 200000.5nm )
				else if (l < surfacePos + 279) {
					previousPoint += 10000.0;
				}
				// Then 20.0um step size (200000.5nm < x < 500000.5)
				else if (l < surfacePos + 294) {
					previousPoint += 20000.0;
				}
				// Then 50.0um step size (500000.5nm < x < 1000000.5)
				else if (l < surfacePos + 304) {
					previousPoint += 50000.0;
				}
				// Then 100.0um step size (1mm < x < 2mm )
				else if (l < surfacePos + 314) {
					previousPoint += 100000.0;
				}
				// Then 200.0um step size (2mm < x < 5mm)
				else if (l < surfacePos + 329) {
					previousPoint += 200000.0;
				}
				// Then 500.0um step size (5mm < x < 10mm)
				else if (l < surfacePos + 339) {
					previousPoint += 500000.0;
				}
				// Then 1.0mm step size (10mm < x)
				else {
					previousPoint += 1000000.0;
				}
			}

			// Get the number of dimensions
			auto dim = opts.getDimensionNumber();
			if (dim > 1) {
				nY = opts.getGridParam(1);
				hY = opts.getGridParam(2);
			}
			if (dim > 2) {
				nZ = opts.getGridParam(3);
				hZ = opts.getGridParam(4);
			}

			return;
		}
		// If it is a geometric gradation grid
		if (gridType == "geometric") {
			// Initialize the value of the previous point
			double previousPoint = 0.0;
			// Set the number of grid points
			nX = opts.getGridParam(0);
			// Set the position of the surface
			auto surfacePos = (IdType)(nX * opts.getVoidPortion() / 100.0);
			// Set the gradation parameters
			double width = 0.1, r = opts.getGridParam(1);

			// Loop on all the grid points
			for (auto l = 0; l <= nX + 1; l++) {
				// Add the previous point
				grid.push_back(previousPoint);
				// 0.1nm step near the surface
				if (l < surfacePos + 1) {
					previousPoint += width;
				}
				else {
					previousPoint += width * pow(r, l - surfacePos - 1);
				}
			}

			// Get the number of dimensions
			auto dim = opts.getDimensionNumber();
			if (dim > 1) {
				nY = opts.getGridParam(2);
				hY = opts.getGridParam(3);
			}
			if (dim > 2) {
				nZ = opts.getGridParam(4);
				hZ = opts.getGridParam(5);
			}

			return;
		}

		throw std::runtime_error("\nThe grid type option was not recognized!");

		return;
	}

	/**
	 * Constructor.
	 *
	 * @param _network The reaction network to use.
	 */
	SolverHandler(NetworkType& _network) :
		network(_network),
		networkName(""),
		nX(0),
		nY(0),
		nZ(0),
		hY(0.0),
		hZ(0.0),
		localXS(0),
		localXM(0),
		localYS(0),
		localYM(0),
		localZS(0),
		localZM(0),
		leftOffset(1),
		rightOffset(1),
		bottomOffset(1),
		topOffset(1),
		frontOffset(1),
		backOffset(1),
		initialVConc(0.0),
		electronicStoppingPower(0.0),
		dimension(-1),
		portion(0.0),
		movingSurface(false),
		bubbleBursting(false),
		isMirror(true),
		useAttenuation(false),
		sputteringYield(0.0),
		fluxHandler(nullptr),
		temperatureHandler(nullptr),
		diffusionHandler(nullptr),
		mutationHandler(nullptr),
		tauBursting(10.0),
		burstingFactor(0.1),
		rngSeed(0),
		heVRatio(4.0),
		previousTime(0.0),
		nXeGB(0.0)
	{
	}

public:
	//! The Constructor
	SolverHandler() = delete;

	~SolverHandler()
	{
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	initializeHandlers(
		std::shared_ptr<core::material::IMaterialHandler> material,
		std::shared_ptr<core::temperature::ITemperatureHandler> tempHandler,
		const options::IOptions& opts) override
	{
		// Determine who I am.
		int myProcId = -1;
		auto xolotlComm = util::getMPIComm();
		MPI_Comm_rank(xolotlComm, &myProcId);

		// Initialize our random number generator.
		bool useRNGSeedFromOptions = false;
		std::tie(useRNGSeedFromOptions, rngSeed) = opts.getRNGSeed();
		if (not useRNGSeedFromOptions) {
			// User didn't give a seed value to use, so
			// use something based on current time and our proc id
			// so that it is different from run to run, and should
			// be different across all processes within a given run.
			rngSeed = time(NULL);
		}
		if (opts.printRNGSeed()) {
			std::cout << "Proc " << myProcId << " using RNG seed value "
					  << rngSeed << std::endl;
		}
		rng = std::unique_ptr<util::RandomNumberGenerator<int, unsigned int>>(
			new util::RandomNumberGenerator<int, unsigned int>(
				rngSeed + myProcId));

		// Set the network loader
		networkName = opts.getNetworkFilename();

		// Set the grid options
		generateGrid(opts);

		// Set the flux handler
		fluxHandler = material->getFluxHandler().get();

		// Set the temperature handler
		temperatureHandler = tempHandler.get();

		// Set the diffusion handler
		diffusionHandler = material->getDiffusionHandler().get();

		// Set the advection handlers
		auto handlers = material->getAdvectionHandler();
		for (auto handler : handlers) {
			advectionHandlers.push_back(handler.get());
		}

		// Set the modified trap-mutation handler
		mutationHandler = material->getTrapMutationHandler().get();

		// Set the minimum size for the average radius computation
		auto numSpecies = network.getSpeciesListSize();
		minRadiusSizes = std::vector<size_t>(numSpecies, 1);
		auto minSizes = opts.getRadiusMinSizes();
		for (auto i = 0; i < std::min(minSizes.size(), minRadiusSizes.size());
			 i++) {
			minRadiusSizes[i] = minSizes[i];
		}

		// Set the initial vacancy concentration
		initialVConc = opts.getInitialVConcentration();

		// Set the electronic stopping power
		electronicStoppingPower = opts.getZeta();

		// Set the number of dimension
		dimension = opts.getDimensionNumber();

		// Set the void portion
		portion = opts.getVoidPortion();

		// Set the sputtering yield
		sputteringYield = opts.getSputteringYield();

		// Set the sputtering yield
		tauBursting = opts.getBurstingDepth();

		// Set the bursting factor
		burstingFactor = opts.getBurstingFactor();

		// Set the HeV ratio
		heVRatio = opts.getHeVRatio();

		// Boundary conditions in the X direction
		if (opts.getMaterial() == "Fuel")
			isMirror = false;

		// Set the boundary conditions (= 1: free surface; = 0: mirror)
		leftOffset = opts.getLeftBoundary();
		rightOffset = opts.getRightBoundary();
		bottomOffset = opts.getBottomBoundary();
		topOffset = opts.getTopBoundary();
		frontOffset = opts.getFrontBoundary();
		backOffset = opts.getBackBoundary();

		// Should we be able to move the surface?
		auto map = opts.getProcesses();
		movingSurface = map["movingSurface"];
		// Should we be able to burst bubbles?
		bubbleBursting = map["bursting"];
		// Should we be able to attenuate the modified trap mutation?
		useAttenuation = map["attenuation"];

		// Some safeguards about what to use with what
		if (leftOffset == 0 &&
			(map["advec"] || map["modifiedTM"] || map["movingSurface"] ||
				map["bursting"])) {
			throw std::runtime_error(
				"\nThe left side of the grid is set to use a reflective "
				"boundary condition but you want to use processes that are "
				"intrinsically related to a free surface (advection, modified "
				"trap mutation, moving surface, bubble bursting).");
		}

		// Complains if processes that should not be used together are used
		if (map["attenuation"] && !map["modifiedTM"]) {
			throw std::runtime_error(
				"\nYou want to use the attenuation on the modified trap "
				"mutation but you are not using the modifiedTM process, it "
				"doesn't make any sense.");
		}
		if (map["modifiedTM"] && !map["reaction"]) {
			throw std::runtime_error(
				"\nYou want to use the modified trap mutation but the reaction "
				"process is not set, it doesn't make any sense.");
		}

		return;
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
	double
	getInitialVConc() const override
	{
		return initialVConc;
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
	double
	getHeVRatio() const override
	{
		return heVRatio;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	createLocalNE(IdType a, IdType b = 1, IdType c = 1) override
	{
		localNE.clear();
		// Create the vector of vectors and fill it with 0.0
		for (auto i = 0; i < a; i++) {
			auto& tempTempVector = localNE.emplace_back();
			for (auto j = 0; j < b; j++) {
				auto& tempVector = tempTempVector.emplace_back();
				for (auto k = 0; k < c; k++) {
					tempVector.push_back({0.0, 0.0, 0.0, 0.0});
				}
			}
		}
	}

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
	burstBubbles() const override
	{
		return bubbleBursting;
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
	core::diffusion::IDiffusionHandler*
	getDiffusionHandler() const override
	{
		return diffusionHandler;
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
	core::modified::ITrapMutationHandler*
	getMutationHandler() const override
	{
		return mutationHandler;
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
	setGBLocation(IdType i, IdType j = 0, IdType k = 0) override
	{
		// Add the coordinates to the GB vector
		if (i >= localXS && i < localXS + std::max(localXM, (IdType)1) &&
			j >= localYS && j < localYS + std::max(localYM, (IdType)1) &&
			k >= localZS && k < localZS + std::max(localZM, (IdType)1)) {
			gbVector.push_back({i, j, k});
		}
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	resetGBVector() override
	{
		gbVector.clear();
	}
};
// end class SolverHandler

} /* namespace handler */
} /* namespace solver */
} /* namespace xolotl */
#endif
