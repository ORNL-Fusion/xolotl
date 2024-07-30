#include <xolotl/factory/perf/PerfHandlerFactory.h>
#include <xolotl/factory/viz/VizHandlerFactory.h>
#include <xolotl/solver/handler/SolverHandler.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/Tokenizer.h>

namespace xolotl
{
namespace solver
{
namespace handler
{
SolverHandler::SolverHandler(NetworkType& _network,
	perf::IPerfHandler& _perfHandler, const options::IOptions& options) :
	network(_network),
	perfHandler(_perfHandler),
	restartFile(options.getRestartFilePath()),
	nX(0),
	nY(0),
	nZ(0),
	hY(0.0),
	hZ(0.0),
	tempGridPower(2.5),
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
	electronicStoppingPower(0.0),
	dimension(-1),
	movingSurface(false),
	bubbleBursting(false),
	isMirror(true),
	isRobin(false),
	useAttenuation(false),
	sameTemperatureGrid(true),
	fluxTempProfile(false),
	sputteringYield(0.0),
	fluxHandler(nullptr),
	temperatureHandler(nullptr),
	vizHandler(factory::viz::VizHandlerFactory::get().generate(options)),
	diffusionHandler(nullptr),
	soretDiffusionHandler(nullptr),
	tauBursting(10.0),
	burstingFactor(0.1),
	rngSeed(0),
	heVRatio(4.0),
	previousTime(0.0),
	nXeGB(0.0),
	gridType(""),
	gridFileName(""),
	gridParam0(0.0),
	gridParam1(0.0),
	gridParam2(0.0),
	gridParam3(0.0),
	gridParam4(0.0),
	gridParam5(0.0)

{
}

SolverHandler::~SolverHandler()
{
}

void
SolverHandler::generateGrid(int surfaceOffset)
{
	// Initial case
	if (surfaceOffset == 0) {
		// Clear the grid
		grid.clear();

		// Doesn't mean anything in 0D
		if (dimension == 0)
			return;

		// Check if we want to read in the grid from a file
		if (gridType == "read") {
			// Open the corresponding file
			std::ifstream inputFile(gridFileName.c_str());
			if (!inputFile)
				std::cerr << "\nCould not open the file containing the grid "
							 "spacing information. "
							 "Aborting!\n"
						  << std::endl;
			// Get the data
			std::string line;
			getline(inputFile, line);

			// Break the line into a vector
			auto tokens = util::Tokenizer<double>{line}();

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

			// Get the number of dimension
			auto dim = dimension;
			if (dim > 1) {
				nY = gridParam0;
				hY = gridParam1;
			}
			if (dim > 2) {
				nZ = gridParam2;
				hZ = gridParam3;
			}

			return;
		}

		// Maybe the user wants a Chebyshev grid
		if (gridType == "cheby") {
			// The first grid point will be at x = 0.0
			grid.push_back(0.0);
			grid.push_back(0.0);

			// Set the number of grid points
			nX = gridParam0;
			auto hx = gridParam1;
			// In that case hx correspond to the full length of the grid
			for (auto l = 1; l <= nX - 1; l++) {
				grid.push_back((hx / 2.0) *
					(1.0 - cos(core::pi * double(l) / double(nX - 1))));
			}
			// The last grid point will be at x = hx
			grid.push_back(hx);

			// Get the number of dimensions
			if (dimension > 1) {
				nY = gridParam2;
				hY = gridParam3;
			}
			if (dimension > 2) {
				nZ = gridParam4;
				hZ = gridParam5;
			}

			return;
		}
		// Check if the user wants a regular grid
		if (gridType == "uniform") {
			// Set the number of grid points
			nX = gridParam0;
			auto hx = gridParam1;
			// The grid will me made of nx + 2 points separated by hx nm
			for (auto l = 0; l <= nX + 1; l++) {
				grid.push_back((double)l * hx);
			}

			// Get the number of dimensions
			if (dimension > 1) {
				nY = gridParam2;
				hY = gridParam3;
			}
			if (dimension > 2) {
				nZ = gridParam4;
				hZ = gridParam5;
			}

			return;
		}

		// If it is not regular do a fine mesh close to the surface and
		// increase the step size when away from the surface
		if (gridType == "nonuniform") {
			// Initialize the value of the previous point
			double previousPoint = 0.0;
			// Set the number of grid points
			nX = gridParam0;
			// Set the position of the surface
			IdType surfacePos = 0;

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
			if (dimension > 1) {
				nY = gridParam1;
				hY = gridParam2;
			}
			if (dimension > 2) {
				nZ = gridParam3;
				hZ = gridParam4;
			}

			return;
		}
		// If it is a geometric gradation grid
		if (gridType == "geometric") {
			// Initialize the value of the previous point
			double previousPoint = 0.0;
			// Set the number of grid points
			nX = gridParam0;
			// Set the gradation parameters
			double width = 0.1, r = gridParam1;

			// Loop on all the grid points
			for (auto l = 0; l <= nX + 1; l++) {
				// Add the previous point
				grid.push_back(previousPoint);
				// 0.1nm step near the surface
				if (l < 1) {
					previousPoint += width;
				}
				else {
					previousPoint += width * pow(r, l - 1);
				}
			}

			// Get the number of dimensions
			if (dimension > 1) {
				nY = gridParam2;
				hY = gridParam3;
			}
			if (dimension > 2) {
				nZ = gridParam4;
				hZ = gridParam5;
			}

			return;
		}

		throw std::runtime_error("\nThe grid type option was not recognized!");
	}
	// Modify grid
	else {
		// Transfer the grid
		oldGrid = grid;

		// Adding grid points case
		if (surfaceOffset > 0) {
			// Compute the distance between what needs to be removed
			double step = grid[surfaceOffset] - grid[0];
			// Check the size of the two used last steps to know if we add a
			// grid point
			double step1 = grid[grid.size() - 2] - grid[grid.size() - 3];
			double step2 = grid[grid.size() - 3] - grid[grid.size() - 4];
			// Modify grid
			if (step1 < step2 - 1.0e-4) {
				grid[grid.size() - 2] += step;
				// Update the last grid point for boundary conditions
				grid[grid.size() - 1] =
					2.0 * grid[grid.size() - 2] - grid[grid.size() - 3];
			}
			// Add the value at the back of the grid
			else {
				// Update the value of the last point
				grid[grid.size() - 1] = grid[grid.size() - 2] + step;
				grid.push_back(
					2.0 * grid[grid.size() - 1] - grid[grid.size() - 2]);
			}
		}
		// Removing grid points case
		// (Work on the current grid because we want to keep the previous
		// geometry in oldGrid)
		else {
			// Compute the distance between what needs to be removed
			double step = grid[-surfaceOffset] - grid[0];
			// Update the value of the last used grid point
			grid[grid.size() - 2] -= step;
			// Remove the last point if it is smaller to the second to last grid
			// point
			while (grid[grid.size() - 2] < grid[grid.size() - 3] + 1.0e-4) {
				grid[grid.size() - 3] = grid[grid.size() - 2];
				grid.pop_back();
			}
			// Update the last grid point for boundary conditions
			grid[grid.size() - 1] =
				2.0 * grid[grid.size() - 2] - grid[grid.size() - 3];
		}
	}

	return;
}

void
SolverHandler::initializeHandlers(core::material::IMaterialHandler* material,
	core::temperature::ITemperatureHandler* tempHandler,
	const options::IOptions& opts)
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
		std::cout << "Proc " << myProcId << " using RNG seed value " << rngSeed
				  << std::endl;
	}
	rng = std::make_unique<util::RandomNumberGenerator<int, unsigned int>>(
		rngSeed + myProcId);

	// Set the flux handler
	fluxHandler = material->getFluxHandler().get();

	// Set the temperature handler
	temperatureHandler = tempHandler;

	// Set the diffusion handlers
	diffusionHandler = material->getDiffusionHandler().get();
	soretDiffusionHandler = material->getSoretDiffusionHandler().get();

	// Set the advection handlers
	auto handlers = material->getAdvectionHandler();
	for (auto handler : handlers) {
		advectionHandlers.push_back(handler.get());
	}

	// Set the minimum size for the average radius computation
	auto numSpecies = network.getSpeciesListSize();
	minRadiusSizes = std::vector<size_t>(numSpecies, 1);
	auto minSizes = opts.getRadiusMinSizes();
	for (auto i = 0; i < std::min(minSizes.size(), minRadiusSizes.size());
		 i++) {
		minRadiusSizes[i] = minSizes[i];
	}

	// Set the initial concentration
	auto initialConcString = opts.getInitialConcentration();
	auto tokens = util::Tokenizer<>{initialConcString}();
	IdType count = 0;
	while (count < tokens.size()) {
		auto comp = std::vector<AmountType>(network.getSpeciesListSize(), 0);

		// Read the cluster type
		auto clusterSpecies = network.parseSpeciesId(tokens[count]);
		// Get the cluster
		comp[clusterSpecies()] = std::stoi(tokens[count + 1]);
		auto clusterId = network.findClusterId(comp);
		// Check that it is present in the network
		if (clusterId == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe requested cluster is not present "
									 "in the network: " +
				tokens[count] + "_" + tokens[count + 1] +
				", cannot use the initial concentration option!");
		}
		else
			initialConc.push_back(std::make_pair<IdType, double>(
				(IdType)clusterId, std::stod(tokens[count + 2])));

		count += 3;
	}

	// Set the electronic stopping power
	electronicStoppingPower = opts.getZeta();

	// Set the number of dimension
	dimension = opts.getDimensionNumber();
	gridType = opts.getGridTypeName();
	gridFileName = opts.getGridFilename();
	gridParam0 = opts.getGridParam(0);
	gridParam1 = opts.getGridParam(1);
	gridParam2 = opts.getGridParam(2);
	gridParam3 = opts.getGridParam(3);
	gridParam4 = opts.getGridParam(4);
	gridParam5 = opts.getGridParam(5);

	// Set the sputtering yield
	sputteringYield = opts.getSputteringYield();

	// Set the sputtering yield
	tauBursting = opts.getBurstingDepth();

	// Set the bursting factor
	burstingFactor = opts.getBurstingFactor();

	// Set the HeV ratio
	heVRatio = opts.getHeVRatio();

	// Which type of temperature grid to use
	if (opts.getTempHandlerName() == "heat" or
		opts.getTempHandlerName() == "ELM") {
		sameTemperatureGrid = false;
		// The temperature grid power only makes sense in this case
		tempGridPower = opts.getTempGridPower();
	}

	// Do we want a flux temporal profile?
	fluxTempProfile = opts.useFluxTimeProfile();

	// Boundary conditions in the X direction
	if (opts.getBCString() == "periodic")
		isMirror = false;
	if (opts.getBCString() == "robin")
		isRobin = true;
	if (opts.getHeatLossPortion() >= 0.0)
		isRobin = true;

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

void
SolverHandler::generateTemperatureGrid()
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
	MPI_Allgatherv(&valuesToSend, toSend, MPI_DOUBLE, &broadcastedTemp, counts,
		displacements, MPI_DOUBLE, xolotlComm);

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
		// Look for it in the temperature grid
		for (auto j = 0; j < nX + 1; j++) {
			double tempLoc1 = 0.0,
				   tempLoc2 = (oldGrid[j] + oldGrid[j + 1]) / 2.0 - oldGrid[1];
			if (j == 0)
				tempLoc1 = oldGrid[0] - oldGrid[1];
			else
				tempLoc1 = (oldGrid[j - 1] + oldGrid[j]) / 2.0 - oldGrid[1];

			if (loc >= tempLoc1 && loc < tempLoc2) {
				double xLoc = (loc - tempLoc1) / (tempLoc2 - tempLoc1);
				double y1 = broadcastedTemp[j], y2 = broadcastedTemp[j + 1];
				toReturn.push_back(y1 + xLoc * (y2 - y1));
				matched = true;
				break;
			}
		}

		if (not matched)
			toReturn.push_back(broadcastedTemp[nX]);
	}

	temperature = toReturn;
}

void
SolverHandler::createLocalNE(IdType a, IdType b, IdType c)
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

bool
SolverHandler::checkForRestart() const
{
    return (not restartFile.empty()) and fs::exists(restartFile);
}

void
SolverHandler::setGBLocation(IdType i, IdType j, IdType k)
{
	// Add the coordinates to the GB vector
	if (i >= localXS && i < localXS + std::max(localXM, (IdType)1) &&
		j >= localYS && j < localYS + std::max(localYM, (IdType)1) &&
		k >= localZS && k < localZS + std::max(localZM, (IdType)1)) {
		gbVector.push_back({i, j, k});
	}
}

std::vector<double>
SolverHandler::interpolateTemperature(std::vector<double> localTemp)
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
	MPI_Allgatherv(&valuesToSend, toSend, MPI_DOUBLE, &broadcastedTemp, counts,
		displacements, MPI_DOUBLE, xolotlComm);

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
		// Look for it in the temperature grid
		for (auto j = 0; j < nX + 1; j++) {
			double tempLoc1 = 0.0,
				   tempLoc2 =
					   (temperatureGrid[j] + temperatureGrid[j + 1]) / 2.0 -
				temperatureGrid[1];
			if (j == 0)
				tempLoc1 = temperatureGrid[0] - temperatureGrid[1];
			else
				tempLoc1 = (temperatureGrid[j - 1] + temperatureGrid[j]) / 2.0 -
					temperatureGrid[1];

			if (loc >= tempLoc1 && loc < tempLoc2) {
				double xLoc = (loc - tempLoc1) / (tempLoc2 - tempLoc1);
				double y1 = broadcastedTemp[j], y2 = broadcastedTemp[j + 1];
				toReturn.push_back(y1 + xLoc * (y2 - y1));
				matched = true;
				break;
			}
		}

		if (not matched) {
			toReturn.push_back(broadcastedTemp[nX + 1]);
		}
	}

	return toReturn;
}
} // namespace handler
} // namespace solver
} // namespace xolotl
