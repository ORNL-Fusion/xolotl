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
	std::vector<std::tuple<int, int, int>> gbVector;

	//! The name of the network file
	std::string networkName;

	//! The name of the free GB file
	std::string gbFileName;

	//! The original network created from the network loader.
	NetworkType& network;

	//! Vector storing the grid in the x direction
	std::vector<double> grid;

	//! The number of grid points in the depth direction.
	int nX;

	//! The number of grid points in the Y direction.
	int nY;

	//! The number of grid points in the Z direction.
	int nZ;

	//! The grid step size in the depth direction.
	double hX;

	//! The grid step size in the y direction.
	double hY;

	//! The grid step size in the z direction.
	double hZ;

	//! The local start of grid points in the X direction.
	int localXS;

	//! The local width of grid points in the X direction.
	int localXM;

	//! The local start of grid points in the Y direction.
	int localYS;

	//! The local width of grid points in the Y direction.
	int localYM;

	//! The local start of grid points in the Z direction.
	int localZS;

	//! The local width of grid points in the Z direction.
	int localZM;

	//! The number of grid points by which the boundary condition should be
	//! shifted at this side.
	int leftOffset, rightOffset, bottomOffset, topOffset, frontOffset,
		backOffset;

	//! The initial vacancy concentration.
	double initialVConc;

	//! The vector of quantities to pass to MOOSE.
	// 0: Xe rate, 1: previous flux, 2: monomer concentration, 3: volume
	// fraction
	std::vector<
		std::vector<std::vector<std::tuple<double, double, double, double>>>>
		localNE;

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

	//! Which type of grid does the used want to use.
	std::string useRegularGrid;

	//! If the user wants to use a Chebyshev grid.
	bool readInGrid;

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

	//! The value to use to seed the random number generator.
	unsigned int rngSeed;

	//! The minimum sizes for average radius computation.
	util::Array<int, 4> minRadiusSizes;

	//! The previous time.
	double previousTime;

	//! The number of xenon atoms that went to the GB
	double nXeGB;

	//! The random number generator to use.
	std::unique_ptr<util::RandomNumberGenerator<int, unsigned int>> rng;

	/**
	 * Method generating the grid in the x direction
	 *
	 * @param nx The number of grid points
	 * @param hx The step size
	 * @param surfacePos The position of the surface on the grid
	 * @param isPSI To know if we want a PSI grid or a NE grid
	 */
	void
	generateGrid(int nx, double hx, int surfacePos)
	{
		// Clear the grid
		grid.clear();

		// Check if we want to read in the grid from a file
		if (readInGrid) {
			// Open the corresponding file
			std::ifstream inputFile(useRegularGrid.c_str());
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
			for (int i = 0; i < tokens.size(); i++) {
				grid.push_back(tokens[i] + offset);
			}

			// Add the last grid point for boundary conditions
			grid.push_back(2.0 * tokens[tokens.size() - 1] -
				tokens[tokens.size() - 2] + offset);

			// Set the number of grid points
			nX = grid.size() - 2;

			return;
		}

		// Maybe the user wants a Chebyshev grid
		if (useRegularGrid == "cheby") {
			// The first grid point will be at x = 0.0
			grid.push_back(0.0);
			grid.push_back(0.0);

			// In that case hx correspond to the full length of the grid
			for (int l = 1; l <= nx - 1; l++) {
				grid.push_back((hx / 2.0) *
					(1.0 - cos(core::pi * double(l) / double(nx - 1))));
			}
			// The last grid point will be at x = hx
			grid.push_back(hx);

			return;
		}
		// Check if the user wants a regular grid
		if (useRegularGrid == "regular") {
			// The grid will me made of nx + 1 points separated by hx nm
			for (int l = 0; l <= nx + 1; l++) {
				grid.push_back((double)l * hx);
			}

			return;
		}
		// If it is not regular do a fine mesh close to the surface and
		// increase the step size when away from the surface
		else if (useRegularGrid == "PSI") {
			// Initialize the value of the previous point
			double previousPoint = 0.0;

			// Loop on all the grid points
			for (int l = 0; l <= nx + 1; l++) {
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

			return;
		}
		// If it is not regular do a fine mesh near points of interests
		else if (useRegularGrid == "NE") {
			// Initialize the value of the previous point
			double previousPoint = 0.0;

			// Loop on all the grid points
			for (int l = 0; l <= nx + 1; l++) {
				// Add the previous point
				grid.push_back(previousPoint);
				// 10nm step near the surface (x < 200nm)
				if (l < surfacePos + 21) {
					previousPoint += 10;
				}
				// 100nm step size (200nm < x < 1um)
				else if (l < surfacePos + 29) {
					previousPoint += 100;
				}
				// 1um step size (1um < x < 5um)
				else if (l < surfacePos + 33) {
					previousPoint += 1000;
				}
				// 5um step size (5um < x < 45um)
				else if (l < surfacePos + 41) {
					previousPoint += 5000;
				}
				// 1um step size (45um < x < 49um)
				else if (l < surfacePos + 45) {
					previousPoint += 1000;
				}
				// 100nm step size
				else if (l < surfacePos + 53) {
					previousPoint += 100;
				}
				// 10nm step size
				else if (l < surfacePos + 93) {
					previousPoint += 10;
				}
				// 100nm step size
				else if (l < surfacePos + 101) {
					previousPoint += 100;
				}
				// 1um step size
				else if (l < surfacePos + 105) {
					previousPoint += 1000;
				}
				// 5um step size
				else if (l < surfacePos + 113) {
					previousPoint += 5000;
				}
				// 1um step size
				else if (l < surfacePos + 117) {
					previousPoint += 1000;
				}
				// 100nm step size
				else if (l < surfacePos + 125) {
					previousPoint += 100;
				}
				// 10nm step size
				else {
					previousPoint += 10;
				}
			}

			return;
		}

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
		hX(0.0),
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
		useRegularGrid(""),
		readInGrid(false),
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
		rngSeed(0),
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
	 * Initialize all the physics handlers that are needed to solve the ADR
	 * equations. \see ISolverHandler.h
	 */
	void
	initializeHandlers(
		std::shared_ptr<core::material::IMaterialHandler> material,
		std::shared_ptr<core::temperature::ITemperatureHandler> tempHandler,
		const options::Options& opts) override
	{
		// Determine who I am.
		int myProcId = -1;
		auto xolotlComm = util::getMPIComm();
		MPI_Comm_rank(xolotlComm, &myProcId);

		// Initialize our random number generator.
		bool useRNGSeedFromOptions = false;
		bool printRNGSeed = false;
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
		// Take the parameter file option by default
		nX = opts.getNX(), nY = opts.getNY(), nZ = opts.getNZ();
		hX = opts.getXStepSize(), hY = opts.getYStepSize(),
		hZ = opts.getZStepSize();
		// Update them if we use an HDF5 file with header group
		if (opts.useHDF5()) {
			int nx = 0, ny = 0, nz = 0;
			double hx = 0.0, hy = 0.0, hz = 0.0;

			io::XFile xfile(networkName);
			auto headerGroup = xfile.getGroup<io::XFile::HeaderGroup>();
			if (headerGroup) {
				headerGroup->read(nx, hx, ny, hy, nz, hz);

				nX = nx, nY = ny, nZ = nz;
				hX = hx, hY = hy, hZ = hz;
			}
		}

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

		// Set the minimum size for the average radius compuation
		minRadiusSizes = opts.getRadiusMinSizes();

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

		// Look at if the user wants to use a regular grid in the x direction
		if (opts.useRegularXGrid())
			useRegularGrid = "regular";
		else if (opts.getMaterial() == "Fuel")
			useRegularGrid = "NE";
		else
			useRegularGrid = "PSI";

		// Boundary conditions in the X direction
		if (opts.getMaterial() == "Fuel")
			isMirror = false;

		// Look at if the user wants to use a Chebyshev grid in the x direction
		if (opts.useChebyshevGrid())
			useRegularGrid = "cheby";

		// Look at if the user wants to read in the grid in the x direction
		if (opts.useReadInGrid()) {
			readInGrid = true;
			useRegularGrid = opts.getGridFilename();
		}

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
			throw std::string("\nThe left side of the grid is set to use a "
							  "reflective boundary condition "
							  "but you want to use processes that are "
							  "intrinsically related to "
							  "a free surface (advection, modified trap "
							  "mutation, moving surface, bubble bursting).");
		}

		// Complains if processes that should not be used together are used
		if (map["attenuation"] && !map["modifiedTM"]) {
			throw std::string("\nYou want to use the attenuation on the "
							  "modified trap mutation "
							  "but you are not using the modifiedTM process, "
							  "it doesn't make any sense.");
		}
		if (map["modifiedTM"] && !map["reaction"]) {
			throw std::string("\nYou want to use the modified trap mutation "
							  "but the reaction process is not set,"
							  " it doesn't make any sense.");
		}

		return;
	}

	/**
	 * Get the grid in the x direction.
	 * \see ISolverHandler.h
	 */
	std::vector<double>
	getXGrid() const override
	{
		return grid;
	}

	/**
	 * Get the step size in the y direction.
	 * \see ISolverHandler.h
	 */
	double
	getStepSizeY() const override
	{
		return hY;
	}

	/**
	 * Get the step size in the z direction.
	 * \see ISolverHandler.h
	 */
	double
	getStepSizeZ() const override
	{
		return hZ;
	}

	/**
	 * Get the number of dimensions of the problem.
	 * \see ISolverHandler.h
	 */
	int
	getDimension() const override
	{
		return dimension;
	}

	/**
	 * Get the initial vacancy concentration.
	 * \see ISolverHandler.h
	 */
	double
	getInitialVConc() const override
	{
		return initialVConc;
	}

	/**
	 * Get the sputtering yield.
	 * \see ISolverHandler.h
	 */
	double
	getSputteringYield() const override
	{
		return sputteringYield;
	}

	/**
	 * Get the depth parameter for bursting.
	 * \see ISolverHandler.h
	 */
	double
	getTauBursting() const override
	{
		return tauBursting;
	}

	/**
	 * Create the local Xe rate vector.
	 * \see ISolverHandler.h
	 */
	void
	createLocalNE(int a, int b = 1, int c = 1) override
	{
		localNE.clear();
		// Create the vector of vectors and fill it with 0.0
		for (int i = 0; i < a; i++) {
			std::vector<std::vector<std::tuple<double, double, double, double>>>
				tempTempVector;
			for (int j = 0; j < b; j++) {
				std::vector<std::tuple<double, double, double, double>>
					tempVector;
				for (int k = 0; k < c; k++) {
					tempVector.push_back(std::make_tuple(0.0, 0.0, 0.0, 0.0));
				}
				tempTempVector.push_back(tempVector);
			}
			localNE.push_back(tempTempVector);
		}
	}

	/**
	 * Set the latest value of the local Xe rate.
	 * \see ISolverHandler.h
	 */
	void
	setLocalXeRate(double rate, int i, int j = 0, int k = 0) override
	{
		std::get<0>(localNE[i][j][k]) += rate;
	}

	/**
	 * Set the whole vector of local Xe rate.
	 * \see ISolverHandler.h
	 */
	void
	setLocalNE(std::vector<
		std::vector<std::vector<std::tuple<double, double, double, double>>>>
			rateVector) override
	{
		localNE = rateVector;
	}

	/**
	 * Get the local Xe rate vector that needs to be passed.
	 * \see ISolverHandler.h
	 */
	std::vector<
		std::vector<std::vector<std::tuple<double, double, double, double>>>>&
	getLocalNE() override
	{
		return localNE;
	}

	/**
	 * Set the latest value of the Xe flux.
	 * \see ISolverHandler.h
	 */
	void
	setPreviousXeFlux(double flux, int i, int j = 0, int k = 0) override
	{
		std::get<1>(localNE[i][j][k]) = flux;
	}

	/**
	 * Set the latest value of the Xe monomer concentration.
	 * \see ISolverHandler.h
	 */
	void
	setMonomerConc(double conc, int i, int j = 0, int k = 0) override
	{
		std::get<2>(localNE[i][j][k]) = conc;
	}

	/**
	 * Set the latest value of the volume fraction.
	 * \see ISolverHandler.h
	 */
	void
	setVolumeFraction(double frac, int i, int j = 0, int k = 0) override
	{
		std::get<3>(localNE[i][j][k]) = frac;
	}

	/**
	 * Set the coordinates covered by the local grid.
	 * \see ISolverHandler.h
	 */
	void
	setLocalCoordinates(
		int xs, int xm, int ys = 0, int ym = 0, int zs = 0, int zm = 0) override
	{
		localXS = xs;
		localXM = xm;
		localYS = ys;
		localYM = ym;
		localZS = zs;
		localZM = zm;

		// Check for free surface boundaries

		// Read the parameter file
		std::ifstream paramFile;
		paramFile.open(gbFileName);
		if (paramFile.good()) {
			// Build an input stream from the string
			util::TokenizedLineReader<int> reader;
			// Get the line
			std::string line;
			getline(paramFile, line);
			auto lineSS = std::make_shared<std::istringstream>(line);
			reader.setInputStream(lineSS);

			// Read the first line
			auto tokens = reader.loadLine();
			// And start looping on the lines
			while (tokens.size() > 0) {
				// Add the coordinates to the GB vector
				setGBLocation(tokens[0], tokens[1], tokens[2]);

				// Read the next line
				getline(paramFile, line);
				lineSS = std::make_shared<std::istringstream>(line);
				reader.setInputStream(lineSS);
				tokens = reader.loadLine();
			}
		}
	}

	/**
	 * Get the coordinates covered by the local grid.
	 * \see ISolverHandler.h
	 */
	void
	getLocalCoordinates(int& xs, int& xm, int& Mx, int& ys, int& ym, int& My,
		int& zs, int& zm, int& Mz) override
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
	 * Get the grid left offset.
	 * \see ISolverHandler.h
	 */
	int
	getLeftOffset() const override
	{
		return leftOffset;
	}

	/**
	 * Get the grid right offset.
	 * \see ISolverHandler.h
	 */
	int
	getRightOffset() const override
	{
		return rightOffset;
	}

	/**
	 * To know if the surface should be able to move.
	 * \see ISolverHandler.h
	 */
	bool
	moveSurface() const override
	{
		return movingSurface;
	}

	/**
	 * To know if the bubble bursting should be used.
	 * \see ISolverHandler.h
	 */
	bool
	burstBubbles() const override
	{
		return bubbleBursting;
	}

	/**
	 * Get the previous time.
	 * \see ISolverHandler.h
	 */
	double
	getPreviousTime() override
	{
		return previousTime;
	}

	/**
	 * Set the previous time.
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
	 * Get the number of Xe that went to the GB.
	 * \see ISolverHandler.h
	 */
	double
	getNXeGB() override
	{
		return nXeGB;
	}

	/**
	 * Set the number of Xe that went to the GB.
	 * \see ISolverHandler.h
	 */
	void
	setNXeGB(double nXe) override
	{
		nXeGB = nXe;
	}

	/**
	 * Get the minimum size for computing average radius.
	 * \see ISolverHandler.h
	 */
	util::Array<int, 4>
	getMinSizes() const override
	{
		return minRadiusSizes;
	}

	/**
	 * Get the flux handler.
	 * \see ISolverHandler.h
	 */
	core::flux::IFluxHandler*
	getFluxHandler() const override
	{
		return fluxHandler;
	}

	/**
	 * Get the temperature handler.
	 * \see ISolverHandler.h
	 */
	core::temperature::ITemperatureHandler*
	getTemperatureHandler() const override
	{
		return temperatureHandler;
	}

	/**
	 * Get the diffusion handler.
	 * \see ISolverHandler.h
	 */
	core::diffusion::IDiffusionHandler*
	getDiffusionHandler() const override
	{
		return diffusionHandler;
	}

	/**
	 * Get the advection handler.
	 * \see ISolverHandler.h
	 */
	core::advection::IAdvectionHandler*
	getAdvectionHandler() const override
	{
		return advectionHandlers[0];
	}

	/**
	 * Get the advection handlers.
	 * \see ISolverHandler.h
	 */
	std::vector<core::advection::IAdvectionHandler*>
	getAdvectionHandlers() const override
	{
		return advectionHandlers;
	}

	/**
	 * Get the modified trap-mutation handler.
	 * \see ISolverHandler.h
	 */
	core::modified::ITrapMutationHandler*
	getMutationHandler() const override
	{
		return mutationHandler;
	}

	/**
	 * Get the network.
	 * \see ISolverHandler.h
	 */
	virtual core::network::IReactionNetwork&
	getNetwork() const override
	{
		return network;
	}

	/**
	 * Get the network name.
	 * \see ISolverHandler.h
	 */
	std::string
	getNetworkName() const override
	{
		return networkName;
	}

	/**
	 * Access the random number generator
	 * The generator will have already been seeded.
	 *
	 * @return The RandomNumberGenerator object to use.
	 */
	util::RandomNumberGenerator<int, unsigned int>&
	getRNG(void) const override
	{
		return *rng;
	}

	/**
	 * Set the file name containing the location of GB.
	 *
	 * @param name The filename
	 */
	void
	setGBFileName(std::string name) override
	{
		gbFileName = name;
	}

	/**
	 * Get the vector containing the location of GB.
	 *
	 * @return The GB vector
	 */
	std::vector<std::tuple<int, int, int>>
	getGBVector() const override
	{
		return gbVector;
	}

	/**
	 * Set the location of one GB grid point.
	 *
	 * @param i, j, k The coordinate of the GB
	 */
	void
	setGBLocation(int i, int j = 0, int k = 0) override
	{
		// Add the coordinates to the GB vector
		if (i >= localXS && i < localXS + std::max(localXM, 1) &&
			j >= localYS && j < localYS + std::max(localYM, 1) &&
			k >= localZS && k < localZS + std::max(localZM, 1)) {
			gbVector.push_back(std::make_tuple(i, j, k));
		}
	}

	/**
	 * Reset the GB vector.
	 */
	void
	resetGBVector() override
	{
		gbVector.clear();
	}

	/**
	 * Get the coordinates covered by the local grid using copying method.
	 * \see ISolverHandler.h
	 */
	void
	getLocalCoordinatesCpy(int* xs, int* xm, int* Mx, int* ys, int* ym, int* My,
		int* zs, int* zm, int* Mz) override
	{
		*xs = localXS;
		*xm = localXM;
		*Mx = nX;
		*ys = localYS;
		*ym = localYM;
		*My = nY;
		*zs = localZS;
		*zm = localZM;
		*Mz = nZ;
	}

	/**
	 * Passing the XeRate at i,j,k point.
	 * \see ISolverHandler.h
	 */
	double
	getXeRatePoint(int i, int j, int k) override
	{
		return std::get<0>(localNE[i][j][k]);
	}

	/**
	 * Passing the XeFlux at i,j,k point.
	 * \see ISolverHandler.h
	 */
	double
	getXeFluxPoint(int i, int j, int k) override
	{
		return std::get<1>(localNE[i][j][k]);
	}

	/**
	 * Passing the XeConc at i,j,k point.
	 * \see ISolverHandler.h
	 */
	double
	getXeConcPoint(int i, int j, int k) override
	{
		return std::get<2>(localNE[i][j][k]);
	}

	/**
	 * Passing the XeVolFrac at i,j,k point.
	 * \see ISolverHandler.h
	 */
	double
	getXeVolFracPoint(int i, int j, int k) override
	{
		return std::get<3>(localNE[i][j][k]);
	}
};
// end class SolverHandler

} /* namespace handler */
} /* namespace solver */
} /* namespace xolotl */
#endif
