#include <petscconf.h>
#include <petscdmda_kokkos.hpp>

#include <xolotl/core/Constants.h>
#include <xolotl/core/Types.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/io/XFile.h>
#include <xolotl/solver/handler/PetscSolver1DHandler.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace solver
{
namespace handler
{
using HostUnmanaged =
	Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

void
PetscSolver1DHandler::createSolverContext(DM& da)
{
	PetscErrorCode ierr;

	// Degrees of freedom is the total number of clusters in the network
	// + moments
	const auto dof = network.getDOF();

	// Set the position of the surface
	surfacePosition = 0;
	if (movingSurface)
		surfacePosition = (IdType)(nX * portion / 100.0);

	// We can update the surface position
	// if we are using a restart file
	if (not networkName.empty() and movingSurface) {
		io::XFile xfile(networkName);
		auto concGroup = xfile.getGroup<io::XFile::ConcentrationGroup>();
		if (concGroup and concGroup->hasTimesteps()) {
			auto tsGroup = concGroup->getLastTimestepGroup();
			assert(tsGroup);
			surfacePosition = tsGroup->readSurface1D();
		}
	}

	// Generate the separate grid for the temperature
	generateTemperatureGrid(surfacePosition);

	// Prints info on one process
	auto xolotlComm = util::getMPIComm();
	int procId = util::getMPIRank();
	if (procId == 0) {
		util::StringStream ss;
		ss << "SolverHandler: 1D simulation with surface BC: ";
		std::string bcString = "periodic";
		if (isMirror)
			bcString = "mirror";
		if (leftOffset == 1)
			ss << "free surface";
		else
			ss << bcString;
		ss << " and bulk BC: ";
		if (rightOffset == 1)
			ss << "free surface";
		else
			ss << bcString;
		for (auto pair : initialConc) {
			ss << ", initial concentration for Id: " << pair.first
			   << " of: " << pair.second << " nm-3";
		}
		ss << ", grid (nm): ";
		for (auto i = 1; i < grid.size() - 1; i++) {
			ss << grid[i] - grid[surfacePosition + 1] << " ";
		}
		ss << std::endl;

		if (not sameTemperatureGrid) {
			ss << "Temperature grid (nm): ";
			for (auto i = 0; i < temperatureGrid.size(); i++) {
				ss << temperatureGrid[i] - temperatureGrid[surfacePosition + 1] << " ";
			}
			ss << std::endl;
		}
		XOLOTL_LOG << ss.str();
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	if (isMirror) {
		ierr = DMDACreate1d(
			xolotlComm, DM_BOUNDARY_MIRROR, nX, dof + 1, 1, NULL, &da);
		checkPetscError(ierr,
			"PetscSolver1DHandler::createSolverContext: "
			"DMDACreate1d failed.");
	}
	else {
		ierr = DMDACreate1d(
			xolotlComm, DM_BOUNDARY_PERIODIC, nX, dof + 1, 1, NULL, &da);
		checkPetscError(ierr,
			"PetscSolver1DHandler::createSolverContext: "
			"DMDACreate1d failed.");
	}
	ierr = DMSetFromOptions(da);
	checkPetscError(ierr,
		"PetscSolver1DHandler::createSolverContext: DMSetFromOptions failed.");
	ierr = DMSetUp(da);
	checkPetscError(
		ierr, "PetscSolver1DHandler::createSolverContext: DMSetUp failed.");
}

void
PetscSolver1DHandler::initializeSolverContext(DM& da, TS& ts)
{
	PetscErrorCode ierr;

	// Degrees of freedom is the total number of clusters in the network
	// + moments
	const auto dof = network.getDOF();

	// Initialize the surface of the first advection handler corresponding to
	// the advection toward the surface (or a dummy one if it is deactivated)
	advectionHandlers[0]->setLocation(grid[surfacePosition + 1] - grid[1]);

	/* The ofill (thought of as a dof by dof 2d (row-oriented) array represents
	 * the nonzero coupling between degrees of freedom at one point with
	 * degrees of freedom on the adjacent point to the left or right.
	 */
	// core::network::IReactionNetwork::SparseFillMap ofill;
	std::vector<core::RowColPair> difEntries;
	std::vector<std::vector<core::RowColPair>> advEntries(
		advectionHandlers.size(), std::vector<core::RowColPair>{});

	// Initialize the temperature handler
	temperatureHandler->initializeTemperature(dof);

	// Fill ofill, the matrix of "off-diagonal" elements that represents
	// diffusion
	diffusionHandler->initialize(network, difEntries);
	// Loop on the advection handlers to account the other "off-diagonal"
	// elements
	for (auto l = 0; l < advectionHandlers.size(); ++l) {
		advectionHandlers[l]->initialize(network, advEntries[l]);
	}

	// Get the local boundaries
	PetscInt xs, xm;
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr,
		"PetscSolver1DHandler::createSolverContext: "
		"DMDAGetCorners failed.");
	// Set it in the handler
	setLocalCoordinates(xs, xm);

	// Tell the network the number of grid points on this process with ghosts
	// TODO: do we need the ghost points?
	network.setGridSize(localXM + 2);

	// Get the diagonal fill
	nNetworkEntries = network.getDiagonalFill(dfill);

	// Load up the block fills
	Mat J;
	ierr = TSGetRHSJacobian(ts, &J, nullptr, nullptr, nullptr);
	checkPetscError(ierr,
		"PetscSolver1DHandler::initializeSolverContext: "
		"TSGetRHSJacobian failed.");
	auto nwEntries = convertToRowColPairList(dof, dfill);
	//
	// "+ 1" for temperature
	auto dSize = localXM * (nNetworkEntries + difEntries.size() + 1);
	// FIXME
	int nAdvec = 0;
	for (auto l = 0; l < advectionHandlers.size(); l++) {
		nAdvec = std::max(nAdvec, advectionHandlers[l]->getNumberOfAdvecting());
	}
	auto oSize =
		localXM * 2 * (difEntries.size() + advEntries.size() * nAdvec + 1);
	auto nPartials = dSize + oSize;
    // TODO: should the count used for this reservation be more exact?
	std::vector<PetscInt> rows, cols;
	rows.reserve(nPartials);
	cols.reserve(nPartials);
	auto mapMatStencilsToCoords =
		[da](const core::RowColPair& component, PetscInt xRow, PetscInt xCol,
			std::vector<PetscInt>& gRows, std::vector<PetscInt>& gCols) {
			PetscFunctionBeginUser;
			MatStencil stCrds[2];
			stCrds[0].i = xRow;
			stCrds[0].c = component[0];
			stCrds[1].i = xCol;
			stCrds[1].c = component[1];
			PetscInt coo[2];
			PetscCall(DMDAMapMatStencilToGlobal(da, 2, stCrds, coo));
			gRows.push_back(coo[0]);
			gCols.push_back(coo[1]);
			PetscFunctionReturn(0);
		};
	std::size_t partialsCount{0};
	for (auto i = localXS; i < localXS + localXM; ++i) {
		// temperature
		mapMatStencilsToCoords({dof, dof}, i, i, rows, cols);
		mapMatStencilsToCoords({dof, dof}, i, i - 1, rows, cols);
		mapMatStencilsToCoords({dof, dof}, i, i + 1, rows, cols);
		partialsCount += 3;
	}
	for (auto i = localXS; i < localXS + localXM; ++i) {
		// diffusion
		for (auto&& component : difEntries) {
			mapMatStencilsToCoords(component, i, i, rows, cols);
			mapMatStencilsToCoords(component, i, i - 1, rows, cols);
			mapMatStencilsToCoords(component, i, i + 1, rows, cols);
		}
		partialsCount += difEntries.size() * 3;
		// advection
		plsm::SpaceVector<double, 3> pos{0.0, 0.0, 0.0};
		pos[0] = (grid[i] + grid[i + 1]) / 2.0 - grid[1];
		for (std::size_t l = 0; l < advectionHandlers.size(); ++l) {
			auto offsets = advectionHandlers[l]->getStencilForAdvection(pos);
			if (advectionHandlers[l]->isPointOnSink(pos)) {
				for (auto&& component : advEntries[l]) {
					mapMatStencilsToCoords(
						component, i, i - offsets[0], rows, cols);
					mapMatStencilsToCoords(
						component, i, i + offsets[0], rows, cols);
				}
			}
			else {
				for (auto&& component : advEntries[l]) {
					mapMatStencilsToCoords(component, i, i, rows, cols);
					mapMatStencilsToCoords(
						component, i, i + offsets[0], rows, cols);
				}
			}
            // Handle potential discrepancies in entry counts
            for (auto i = advEntries[l].size(); i < nAdvec; ++i) {
                rows.push_back(-1);
                cols.push_back(-1);
            }
			partialsCount += nAdvec * 2;
		}
		// network
		for (auto&& component : nwEntries) {
			mapMatStencilsToCoords(component, i, i, rows, cols);
		}
		partialsCount += nwEntries.size();
	}
	//
	std::cout << "count: " << partialsCount << "\nnPartialsEst: " << nPartials
			  << "\nnPartials: " << rows.size() << std::endl;
    nPartials = rows.size();
	ierr = MatSetPreallocationCOO(J, nPartials, rows.data(), cols.data());
	checkPetscError(ierr,
		"PetscSolver1DHandler::initializeSolverContext: "
		"MatSetPreallocationCOO failed.");

	// Initialize the arrays for the reaction partial derivatives
	vals = Kokkos::View<double*>("solverPartials", nPartials);

	// Set the size of the partial derivatives vectors
	reactingPartialsForCluster.resize(dof, 0.0);

	// Initialize the flux handler
	fluxHandler->initializeFluxHandler(network, surfacePosition, grid);
}

void
PetscSolver1DHandler::initializeConcentration(DM& da, Vec& C)
{
	PetscErrorCode ierr;

	// Pointer for the concentration vector
	PetscScalar** concentrations = nullptr;
	ierr = DMDAVecGetArrayDOF(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver1DHandler::initializeConcentration: "
		"DMDAVecGetArrayDOF failed.");

	// Initialize the last temperature at each grid point on this process
	for (auto i = 0; i < localXM + 2; i++) {
		temperature.push_back(0.0);
	}

	// Get the last time step written in the HDF5 file
	bool hasConcentrations = false;
	std::unique_ptr<io::XFile> xfile;
	std::unique_ptr<io::XFile::ConcentrationGroup> concGroup;
	if (not networkName.empty()) {
		xfile = std::make_unique<io::XFile>(networkName);
		concGroup = xfile->getGroup<io::XFile::ConcentrationGroup>();
		hasConcentrations = (concGroup and concGroup->hasTimesteps());
	}

	// Give the surface position to the temperature handler
	temperatureHandler->updateSurfacePosition(surfacePosition);

	// Initialize the grid for the diffusion
	diffusionHandler->initializeDiffusionGrid(
		advectionHandlers, grid, localXM, localXS);

	// Initialize the grid for the advection
	advectionHandlers[0]->initializeAdvectionGrid(
		advectionHandlers, grid, localXM, localXS);

	// Pointer for the concentration vector at a specific grid point
	PetscScalar* concOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	// + moments
	const auto dof = network.getDOF();

	// Loop on all the grid points
	for (auto i = (PetscInt)localXS - 1;
		 i <= (PetscInt)localXS + (PetscInt)localXM; i++) {
		// Temperature
		plsm::SpaceVector<double, 3> gridPosition{0.0, 0.0, 0.0};
		if (i < 0)
			gridPosition[0] =
				(temperatureGrid[0] - temperatureGrid[surfacePosition + 1]) /
				(temperatureGrid[temperatureGrid.size() - 1] -
					temperatureGrid[surfacePosition + 1]);
		else
			gridPosition[0] =
				((temperatureGrid[i] + temperatureGrid[i + 1]) / 2.0 -
					temperatureGrid[surfacePosition + 1]) /
				(temperatureGrid[temperatureGrid.size() - 1] -
					temperatureGrid[surfacePosition + 1]);
		auto temp = temperatureHandler->getTemperature(gridPosition, 0.0);
		temperature[i - localXS + 1] = temp;

		// Boundary conditions
		if (i < localXS || i >= localXS + localXM)
			continue;

		concOffset = concentrations[i];
		concOffset[dof] = temp;

		// Loop on all the clusters to initialize at 0.0
		for (auto n = 0; n < dof; n++) {
			concOffset[n] = 0.0;
		}

		// Initialize the option specified concentration
		if (i >= surfacePosition + leftOffset and not hasConcentrations and
			i < nX - rightOffset) {
			for (auto pair : initialConc) {
				concOffset[pair.first] = pair.second;
			}
		}
	}

	// If the concentration must be set from the HDF5 file
	if (hasConcentrations) {
		// Read the concentrations from the HDF5 file for
		// each of our grid points.
		assert(concGroup);
		auto tsGroup = concGroup->getLastTimestepGroup();
		assert(tsGroup);
		auto myConcs = tsGroup->readConcentrations(*xfile, localXS, localXM);

		// Apply the concentrations we just read.
		for (auto i = 0; i < localXM; ++i) {
			concOffset = concentrations[localXS + i];

			for (auto const& currConcData : myConcs[i]) {
				concOffset[currConcData.first] = currConcData.second;
			}
			// Get the temperature
			double temp = myConcs[i][myConcs[i].size() - 1].second;
			temperature[i + 1] = temp;
		}
	}

	// Update the network with the temperature
	auto networkTemp = interpolateTemperature(surfacePosition);
	std::vector<double> depths;
	for (auto i = 0; i < networkTemp.size(); i++) {
		if (localXS + i == nX + 1)
			depths.push_back(grid[localXS + i] - grid[surfacePosition + 1]);
		else
			depths.push_back((grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
				grid[surfacePosition + 1]);
	}
	network.setTemperatures(networkTemp, depths);

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver1DHandler::initializeConcentration: "
		"DMDAVecRestoreArrayDOF failed.");

	return;
}

void
PetscSolver1DHandler::initGBLocation(DM& da, Vec& C)
{
	PetscErrorCode ierr;

	// Pointer for the concentration vector
	PetscScalar** concentrations = nullptr;
	ierr = DMDAVecGetArrayDOF(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver1DHandler::initGBLocation: "
		"DMDAVecGetArrayDOF failed.");

	// Pointer for the concentration vector at a specific grid point
	PetscScalar* concOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	// + moments
	const auto dof = network.getDOF();

	// Need to use the NE network here
	using NetworkType = core::network::NEReactionNetwork;
	using Spec = typename NetworkType::Species;
	auto& neNetwork = dynamic_cast<NetworkType&>(network);

	// Loop on the GB
	for (auto const& pair : gbVector) {
		// Get the coordinate of the point
		auto xi = std::get<0>(pair);
		// Check if we are on the right process
		if (xi >= localXS && xi < localXS + localXM) {
			// Get the local concentration
			concOffset = concentrations[xi];

			auto hConcs = HostUnmanaged(concOffset, dof);
			auto dConcs = Kokkos::View<double*>("Concentrations", dof);
			deep_copy(dConcs, hConcs);

			// Transfer the local amount of Xe clusters
			setLocalXeRate(
				neNetwork.getTotalAtomConcentration(dConcs, Spec::Xe, 1),
				xi - localXS);

			// Loop on all the clusters to initialize at 0.0
			for (auto n = 0; n < dof; n++) {
				concOffset[n] = 0.0;
			}
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver1DHandler::initGBLocation: "
		"DMDAVecRestoreArrayDOF failed.");

	return;
}

std::vector<std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
PetscSolver1DHandler::getConcVector(DM& da, Vec& C)
{
	// Initial declaration
	PetscErrorCode ierr;
	const double* gridPointSolution = nullptr;

	// Pointer for the concentration vector
	PetscScalar** concentrations = nullptr;
	ierr = DMDAVecGetArrayDOFRead(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver1DHandler::getConcVector: "
		"DMDAVecGetArrayDOFRead failed.");

	// Get the network and dof
	auto& network = getNetwork();
	const auto dof = network.getDOF();

	// Create the vector for the concentrations
	std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
		toReturn;
	std::vector<std::vector<std::pair<IdType, double>>> tempTempVector;

	// Loop on the grid points
	for (auto i = 0; i < localXM; ++i) {
		gridPointSolution = concentrations[localXS + i];

		// Create the temporary vector for this grid point
		std::vector<std::pair<IdType, double>> tempVector;
		for (auto l = 0; l < dof + 1; ++l) {
			if (std::fabs(gridPointSolution[l]) > 1.0e-16) {
				tempVector.push_back(std::make_pair(l, gridPointSolution[l]));
			}
		}
		tempTempVector.push_back(tempVector);
	}
	std::vector<std::vector<std::vector<std::pair<IdType, double>>>>
		tempTempTempVector;
	tempTempTempVector.push_back(tempTempVector);
	toReturn.push_back(tempTempTempVector);

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver1DHandler::getConcVector: "
		"DMDAVecRestoreArrayDOFRead failed.");

	return toReturn;
}

void
PetscSolver1DHandler::setConcVector(DM& da, Vec& C,
	std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>&
		concVector)
{
	PetscErrorCode ierr;

	// Pointer for the concentration vector
	PetscScalar* gridPointSolution = nullptr;
	PetscScalar** concentrations = nullptr;
	ierr = DMDAVecGetArrayDOF(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver1DHandler::setConcVector: "
		"DMDAVecGetArrayDOF failed.");

	// Loop on the grid points
	for (auto i = 0; i < localXM; ++i) {
		gridPointSolution = concentrations[localXS + i];

		// Loop on the given vector
		for (auto l = 0; l < concVector[0][0][i].size(); l++) {
			gridPointSolution[concVector[0][0][i][l].first] =
				concVector[0][0][i][l].second;
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver1DHandler::setConcVector: "
		"DMDAVecRestoreArrayDOF failed.");

	// Get the complete data array, including ghost cells to set the temperature
	// at the ghost points
	Vec localSolution;
	ierr = DMGetLocalVector(da, &localSolution);
	checkPetscError(ierr,
		"PetscSolver1DHandler::setConcVector: "
		"DMGetLocalVector failed.");
	ierr = DMGlobalToLocalBegin(da, C, INSERT_VALUES, localSolution);
	checkPetscError(ierr,
		"PetscSolver1DHandler::setConcVector: "
		"DMGlobalToLocalBegin failed.");
	ierr = DMGlobalToLocalEnd(da, C, INSERT_VALUES, localSolution);
	checkPetscError(ierr,
		"PetscSolver1DHandler::setConcVector: "
		"DMGlobalToLocalEnd failed.");
	// Get the array of concentration
	ierr = DMDAVecGetArrayDOFRead(da, localSolution, &concentrations);
	checkPetscError(ierr,
		"PetscSolver1DHandler::setConcVector: "
		"DMDAVecGetArrayDOFRead failed.");

	// Getthe DOF of the network
	const auto dof = network.getDOF();

	// Loop on the grid points
	for (auto i = -1; i <= localXM; ++i) {
		gridPointSolution = concentrations[localXS + i];

		// Get the temperature
		temperature[i + 1] = gridPointSolution[dof];
	}
	// Update the network with the temperature
	auto networkTemp = interpolateTemperature(surfacePosition);
	std::vector<double> depths;
	for (auto i = 0; i < networkTemp.size(); i++) {
		if (localXS + i == nX + 1)
			depths.push_back(grid[localXS + i] - grid[surfacePosition + 1]);
		else
			depths.push_back((grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
				grid[surfacePosition + 1]);
	}
	network.setTemperatures(networkTemp, depths);

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, localSolution, &concentrations);
	checkPetscError(ierr,
		"PetscSolver1DHandler::setConcVector: "
		"DMDAVecRestoreArrayDOFRead failed.");
	ierr = DMRestoreLocalVector(da, &localSolution);
	checkPetscError(ierr,
		"PetscSolver1DHandler::setConcVector: "
		"DMRestoreLocalVector failed.");

	return;
}

void
PetscSolver1DHandler::updateConcentration(
	TS& ts, Vec& localC, Vec& F, PetscReal ftime)
{
	PetscErrorCode ierr;

	// Get the local data vector from PETSc
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr,
		"PetscSolver1DHandler::updateConcentration: "
		"TSGetDM failed.");

	// Pointers to the PETSc arrays that start at the beginning (localXS) of the
	// local array!
	PetscOffsetView<const PetscScalar**> concs;
	ierr = DMDAVecGetKokkosOffsetViewDOF(da, localC, &concs);
	checkPetscError(ierr,
		"PetscSolver1DHandler::updateConcentration: "
		"DMDAVecGetKokkosOffsetViewDOF (localC) failed.");
	PetscOffsetView<PetscScalar**> updatedConcs;
	ierr = DMDAVecGetKokkosOffsetViewDOFWrite(da, F, &updatedConcs);
	checkPetscError(ierr,
		"PetscSolver1DHandler::updateConcentration: "
		"DMDAVecGetKokkosOffsetViewDOFWrite (F) failed.");

	// Degrees of freedom is the total number of clusters in the network
	const auto dof = network.getDOF();

	// Computing the trapped atom concentration is only needed for the
	// attenuation
	if (useAttenuation) {
		// Compute the total concentration of atoms contained in bubbles
		double atomConc = 0.0;

		auto& psiNetwork =
			dynamic_cast<core::network::IPSIReactionNetwork&>(network);

		// Loop over grid points to get the atom concentration
		// near the surface
		for (auto xi = localXS; xi < localXS + localXM; xi++) {
			// Boundary conditions
			if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset)
				continue;

			// We are only interested in the helium near the surface
			if ((grid[xi] + grid[xi + 1]) / 2.0 - grid[surfacePosition + 1] >
				2.0)
				continue;

			// Get the concentrations at this grid point
			auto concOffset = subview(concs, xi, Kokkos::ALL).view();

			// Sum the total atom concentration
			atomConc +=
				psiNetwork.getTotalTrappedHeliumConcentration(concOffset, 0) *
				(grid[xi + 1] - grid[xi]);
		}

		// Share the concentration with all the processes
		double totalAtomConc = 0.0;
		auto xolotlComm = util::getMPIComm();
		MPI_Allreduce(
			&atomConc, &totalAtomConc, 1, MPI_DOUBLE, MPI_SUM, xolotlComm);

		// Set the disappearing rate in the modified TM handler
		psiNetwork.updateTrapMutationDisappearingRate(totalAtomConc);
	}

	// Declarations for variables used in the loop
	using ConcSubView = Kokkos::View<const double*>;
	Kokkos::Array<ConcSubView, 3> concVector;
	plsm::SpaceVector<double, 3> gridPosition{0.0, 0.0, 0.0};

	// Loop over grid points first for the temperature, including the ghost
	// points
	bool tempHasChanged = false;
	for (auto xi = (PetscInt)localXS - 1;
		 xi <= (PetscInt)localXS + (PetscInt)localXM; xi++) {
		// Heat condition
		if (xi == surfacePosition && xi >= localXS && xi < localXS + localXM) {
			// Compute the old and new array offsets
			auto concOffset = subview(concs, xi, Kokkos::ALL).view();
			auto updatedConcOffset =
				subview(updatedConcs, xi, Kokkos::ALL).view();

			// Fill the concVector with the pointer to the middle, left, and
			// right grid points
			concVector[0] = concOffset; // middle
			concVector[1] =
				subview(concs, (PetscInt)xi - 1, Kokkos::ALL).view(); // left
			concVector[2] = subview(concs, xi + 1, Kokkos::ALL).view(); // right

			// Compute the left and right hx
			double hxLeft = 0.0, hxRight = 0.0;
			if (xi >= 1 && xi < nX) {
				hxLeft =
					(temperatureGrid[xi + 1] - temperatureGrid[xi - 1]) / 2.0;
				hxRight = (temperatureGrid[xi + 2] - temperatureGrid[xi]) / 2.0;
			}
			else if (xi < 1) {
				hxLeft = temperatureGrid[xi + 1] - temperatureGrid[xi];
				hxRight = (temperatureGrid[xi + 2] - temperatureGrid[xi]) / 2.0;
			}
			else {
				hxLeft =
					(temperatureGrid[xi + 1] - temperatureGrid[xi - 1]) / 2.0;
				hxRight = temperatureGrid[xi + 1] - temperatureGrid[xi];
			}

			temperatureHandler->computeTemperature(
				concVector.data(), updatedConcOffset, hxLeft, hxRight, xi);
		}

		// Compute the old and new array offsets
		auto concOffset = subview(concs, xi, Kokkos::ALL).view();
		auto updatedConcOffset = subview(updatedConcs, xi, Kokkos::ALL).view();

		// Set the grid fraction
		if (xi < 0) {
			gridPosition[0] = (grid[0] - grid[surfacePosition + 1]) /
				(grid.back() - grid[surfacePosition + 1]);
		}
		else {
			gridPosition[0] =
				((grid[xi] + grid[xi + 1]) / 2.0 - grid[surfacePosition + 1]) /
				(grid.back() - grid[surfacePosition + 1]);
		}

		// Get the temperature from the temperature handler
		temperatureHandler->setTemperature(concOffset);
		double temp = temperatureHandler->getTemperature(gridPosition, ftime);

		// Update the network if the temperature changed
		if (std::fabs(temperature[xi + 1 - localXS] - temp) > 0.1) {
			temperature[xi + 1 - localXS] = temp;
			tempHasChanged = true;
		}

		// Boundary conditions
		// Everything to the left of the surface is empty
		if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset) {
			continue;
		}
		// Free surface GB
		bool skip = false;
		for (auto& pair : gbVector) {
			if (xi == std::get<0>(pair)) {
				skip = true;
				break;
			}
		}
		if (skip)
			continue;

		// Fill the concVector with the pointer to the middle, left, and right
		// grid points
		concVector[0] = concOffset; // middle
		concVector[1] =
			subview(concs, (PetscInt)xi - 1, Kokkos::ALL).view(); // left
		concVector[2] = subview(concs, xi + 1, Kokkos::ALL).view(); // right

		// Compute the left and right hx
		double hxLeft = 0.0, hxRight = 0.0;
		if (xi >= 1 && xi < nX) {
			hxLeft = (temperatureGrid[xi + 1] - temperatureGrid[xi - 1]) / 2.0;
			hxRight = (temperatureGrid[xi + 2] - temperatureGrid[xi]) / 2.0;
		}
		else if (xi < 1) {
			hxLeft = temperatureGrid[xi + 1] - temperatureGrid[xi];
			hxRight = (temperatureGrid[xi + 2] - temperatureGrid[xi]) / 2.0;
		}
		else {
			hxLeft = (temperatureGrid[xi + 1] - temperatureGrid[xi - 1]) / 2.0;
			hxRight = temperatureGrid[xi + 1] - temperatureGrid[xi];
		}

		// ---- Compute the temperature over the locally owned part of the grid
		// -----
		if (xi >= localXS && xi < localXS + localXM) {
			temperatureHandler->computeTemperature(
				concVector.data(), updatedConcOffset, hxLeft, hxRight, xi);
		}
	}

	// Share the information with all the processes
	bool totalTempHasChanged = false;
	auto xolotlComm = util::getMPIComm();
	MPI_Allreduce(&tempHasChanged, &totalTempHasChanged, 1, MPI_C_BOOL, MPI_LOR,
		xolotlComm);

	if (totalTempHasChanged) {
		// Update the network with the temperature
		auto networkTemp = interpolateTemperature(surfacePosition);
		std::vector<double> depths;
		for (auto i = 0; i < networkTemp.size(); i++) {
			if (localXS + i == nX + 1)
				depths.push_back(grid[localXS + i] - grid[surfacePosition + 1]);
			else
				depths.push_back(
					(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
					grid[surfacePosition + 1]);
		}
		network.setTemperatures(networkTemp, depths);
	}

	// Loop over grid points computing ODE terms for each grid point
	for (auto xi = localXS; xi < localXS + localXM; xi++) {
		// Compute the old and new array offsets
		auto concOffset = subview(concs, xi, Kokkos::ALL).view();
		auto updatedConcOffset = subview(updatedConcs, xi, Kokkos::ALL).view();

		// Fill the concVector with the pointer to the middle, left, and right
		// grid points
		concVector[0] = concOffset; // middle
		concVector[1] =
			subview(concs, (PetscInt)xi - 1, Kokkos::ALL).view(); // left
		concVector[2] = subview(concs, xi + 1, Kokkos::ALL).view(); // right

		// Compute the left and right hx
		double hxLeft = 0.0, hxRight = 0.0;
		if (xi >= 1 && xi < nX) {
			hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
			hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
		}
		else if (xi < 1) {
			hxLeft = grid[xi + 1] - grid[xi];
			hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
		}
		else {
			hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
			hxRight = grid[xi + 1] - grid[xi];
		}

		// Boundary conditions
		// Everything to the left of the surface is empty
		if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset) {
			continue;
		}
		// Free surface GB
		bool skip = false;
		for (auto& pair : gbVector) {
			if (xi == std::get<0>(pair)) {
				skip = true;
				break;
			}
		}
		if (skip)
			continue;

		// ----- Account for flux of incoming particles -----
		fluxHandler->computeIncidentFlux(
			ftime, updatedConcOffset, xi, surfacePosition);

		// ---- Compute diffusion over the locally owned part of the grid -----
		diffusionHandler->computeDiffusion(network,
			core::StencilConcArray{concVector.data(), 3}, updatedConcOffset,
			hxLeft, hxRight, xi - localXS);

		// ---- Compute advection over the locally owned part of the grid -----
		// Set the grid position
		gridPosition[0] = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];
		for (auto i = 0; i < advectionHandlers.size(); i++) {
			advectionHandlers[i]->computeAdvection(network, gridPosition,
				core::StencilConcArray{concVector.data(), 3}, updatedConcOffset,
				hxLeft, hxRight, xi - localXS);
		}

		auto surfacePos = grid[surfacePosition + 1];
		auto curXPos = (grid[xi] + grid[xi + 1]) / 2.0;
		auto prevXPos = (grid[xi - 1] + grid[xi]) / 2.0;
		auto curDepth = curXPos - surfacePos;
		auto curSpacing = curXPos - prevXPos;

		// ----- Compute the reaction fluxes over the locally owned part of the
		// grid -----
		fluxCounter->increment();
		fluxTimer->start();
		network.computeAllFluxes(concOffset, updatedConcOffset,
			xi + 1 - localXS, curDepth, curSpacing);
		fluxTimer->stop();
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreKokkosOffsetViewDOF(da, localC, &concs);
	checkPetscError(ierr,
		"PetscSolver1DHandler::updateConcentration: "
		"DMDAVecRestoreKokkosOffsetViewDOF (localC) failed.");
	ierr = DMDAVecRestoreKokkosOffsetViewDOFWrite(da, F, &updatedConcs);
	checkPetscError(ierr,
		"PetscSolver1DHandler::updateConcentration: "
		"DMDAVecRestoreKokkosOffsetViewDOFWrite (F) failed.");
}

void
PetscSolver1DHandler::computeJacobian(
	TS& ts, Vec& localC, Mat& J, PetscReal ftime)
{
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr,
		"PetscSolver1DHandler::computeJacobian: "
		"TSGetDM failed.");

	PetscOffsetView<const PetscScalar**> concs;
	ierr = DMDAVecGetKokkosOffsetViewDOF(da, localC, &concs);
	checkPetscError(ierr,
		"PetscSolver1DHandler::computeJacobian: "
		"DMDAVecGetKokkosOffsetViewDOF failed.");

	// Degrees of freedom is the total number of clusters in the network
	const auto dof = network.getDOF();

	// Get the total number of diffusing clusters
	const auto nDiff = std::max(diffusionHandler->getNumberOfDiffusing(), 0);

	// Get the total number of advecting clusters
	int nAdvec = 0;
	for (auto l = 0; l < advectionHandlers.size(); l++) {
		nAdvec = std::max(nAdvec, advectionHandlers[l]->getNumberOfAdvecting());
	}

	// Arguments for MatSetValuesStencil called below
	MatStencil row, cols[3];
	PetscScalar tempVals[3];
	IdType tempIndices[1];
	auto diffVals = std::vector<PetscScalar>(3 * nDiff);
	auto diffIndices = std::vector<IdType>(nDiff);
	auto advecVals = std::vector<PetscScalar>(2 * nAdvec);
	auto advecIndices = std::vector<IdType>(nAdvec);
	plsm::SpaceVector<double, 3> gridPosition{0.0, 0.0, 0.0};

	/*
	 Loop over grid points for the temperature, including ghosts
	 */
	bool tempHasChanged = false;
	auto hTempVals = Kokkos::View<double*, Kokkos::HostSpace>(
		"Host Temp Jac Vals", localXM * 3);
	std::size_t valIndex = 0;
	for (auto xi = (PetscInt)localXS - 1;
		 xi <= (PetscInt)localXS + (PetscInt)localXM; xi++) {
		// Compute the left and right hx
		double hxLeft = 0.0, hxRight = 0.0;
		if (xi >= 1 && xi < nX) {
			hxLeft = (temperatureGrid[xi + 1] - temperatureGrid[xi - 1]) / 2.0;
			hxRight = (temperatureGrid[xi + 2] - temperatureGrid[xi]) / 2.0;
		}
		else if (xi == -1) {
			hxLeft = 0.0;
			hxRight = (temperatureGrid[xi + 2] + temperatureGrid[xi + 1]) / 2.0;
		}
		else if (xi < 1) {
			hxLeft = temperatureGrid[xi + 1] - temperatureGrid[xi];
			hxRight = (temperatureGrid[xi + 2] - temperatureGrid[xi]) / 2.0;
		}
		else {
			hxLeft = (temperatureGrid[xi + 1] - temperatureGrid[xi - 1]) / 2.0;
			hxRight = temperatureGrid[xi + 1] - temperatureGrid[xi];
		}

		auto tempIndex = valIndex;
		if (xi >= localXS && xi < localXS + localXM) {
			valIndex += 3;
		}

		// Heat condition
		if (xi == surfacePosition && xi >= localXS && xi < localXS + localXM) {
			// Get the partial derivatives for the temperature
			auto setValues = temperatureHandler->computePartialsForTemperature(
				tempVals, tempIndices, hxLeft, hxRight, xi);

			if (setValues) {
				hTempVals(tempIndex + 0) += tempVals[0];
				hTempVals(tempIndex + 1) += tempVals[1];
				hTempVals(tempIndex + 2) += tempVals[2];
			}
		}

		// Get the concentrations at this grid point
		auto concOffset = subview(concs, xi, Kokkos::ALL).view();

		// Set the grid fraction
		if (xi < 0) {
			gridPosition[0] =
				(temperatureGrid[0] - temperatureGrid[surfacePosition + 1]) /
				(temperatureGrid.back() - temperatureGrid[surfacePosition + 1]);
		}
		else {
			gridPosition[0] =
				((temperatureGrid[xi] + temperatureGrid[xi + 1]) / 2.0 -
					temperatureGrid[surfacePosition + 1]) /
				(temperatureGrid.back() - temperatureGrid[surfacePosition + 1]);
		}

		// Get the temperature from the temperature handler
		temperatureHandler->setTemperature(concOffset);
		double temp = temperatureHandler->getTemperature(gridPosition, ftime);

		// Update the network if the temperature changed
		if (std::fabs(temperature[xi + 1 - localXS] - temp) > 0.1) {
			temperature[xi + 1 - localXS] = temp;
			tempHasChanged = true;
		}

		// Boundary conditions
		// Everything to the left of the surface is empty
		if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset)
			continue;
		// Free surface GB
		bool skip = false;
		for (auto& pair : gbVector) {
			if (xi == std::get<0>(pair)) {
				skip = true;
				break;
			}
		}
		if (skip)
			continue;

		// Get the partial derivatives for the temperature
		if (xi >= localXS && xi < localXS + localXM) {
			auto setValues = temperatureHandler->computePartialsForTemperature(
				tempVals, tempIndices, hxLeft, hxRight, xi);

			if (setValues) {
				hTempVals(tempIndex + 0) += tempVals[0];
				hTempVals(tempIndex + 1) += tempVals[1];
				hTempVals(tempIndex + 2) += tempVals[2];
			}
		}
	}
	deep_copy(subview(vals, std::make_pair(IdType{0}, localXM * 3)), hTempVals);

	// Share the information with all the processes
	bool totalTempHasChanged = false;
	auto xolotlComm = util::getMPIComm();
	MPI_Allreduce(&tempHasChanged, &totalTempHasChanged, 1, MPI_C_BOOL, MPI_LOR,
		xolotlComm);

	if (totalTempHasChanged) {
		// Update the network with the temperature
		auto networkTemp = interpolateTemperature(surfacePosition);
		std::vector<double> depths;
		for (auto i = 0; i < networkTemp.size(); i++) {
			if (localXS + i == nX + 1)
				depths.push_back(grid[localXS + i] - grid[surfacePosition + 1]);
			else
				depths.push_back(
					(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
					grid[surfacePosition + 1]);
		}
		network.setTemperatures(networkTemp, depths);
	}

	// Computing the trapped atom concentration is only needed for the
	// attenuation
	if (useAttenuation) {
		// Compute the total concentration of atoms contained in bubbles
		double atomConc = 0.0;

		auto& psiNetwork =
			dynamic_cast<core::network::IPSIReactionNetwork&>(network);

		// Loop over grid points to get the atom concentration
		// near the surface
		for (auto xi = localXS; xi < localXS + localXM; xi++) {
			// Boundary conditions
			if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset)
				continue;

			// We are only interested in the helium near the surface
			if ((grid[xi] + grid[xi + 1]) / 2.0 - grid[surfacePosition + 1] >
				2.0)
				continue;

			// Get the concentrations at this grid point
			auto concOffset = subview(concs, xi, Kokkos::ALL).view();

			// Sum the total atom concentration
			atomConc +=
				psiNetwork.getTotalTrappedHeliumConcentration(concOffset, 0) *
				(grid[xi + 1] - grid[xi]);
		}

		// Share the concentration with all the processes
		double totalAtomConc = 0.0;
		auto xolotlComm = util::getMPIComm();
		MPI_Allreduce(
			&atomConc, &totalAtomConc, 1, MPI_DOUBLE, MPI_SUM, xolotlComm);

		// Set the disappearing rate in the modified TM handler
		psiNetwork.updateTrapMutationDisappearingRate(totalAtomConc);
	}

	// Loop over the grid points
	for (auto xi = localXS; xi < localXS + localXM; xi++) {
		// Boundary conditions
		// Everything to the left of the surface is empty
		if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset) {
            valIndex += 3 * nDiff;
			valIndex += 2 * nAdvec;
            valIndex += nNetworkEntries;
			continue;
        }

		// Free surface GB
        if (std::find_if(begin(gbVector), end(gbVector), [=](auto&& pair) {
                    return xi == pair[0]; }) != end(gbVector)) {
            // TODO: If the gbVector is initialized before the preallocation, we
            // could simply avoid the extra entries
            valIndex += 3 * nDiff;
			valIndex += 2 * nAdvec;
            valIndex += nNetworkEntries;
            continue;
        }

		// Compute the left and right hx
		double hxLeft = 0.0, hxRight = 0.0;
		if (xi >= 1 && xi < nX) {
			hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
			hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
		}
		else if (xi < 1) {
			hxLeft = grid[xi + 1] - grid[xi];
			hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
		}
		else {
			hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
			hxRight = grid[xi + 1] - grid[xi];
		}

		// Get the partial derivatives for the diffusion
		diffusionHandler->computePartialsForDiffusion(network, diffVals.data(),
			diffIndices.data(), hxLeft, hxRight, xi - localXS);

		auto hDiffVals = HostUnmanaged(diffVals.data(), 3 * nDiff);
		deep_copy(subview(vals, std::make_pair(valIndex, valIndex + 3 * nDiff)),
			hDiffVals);
		valIndex += 3 * nDiff;

		// Get the partial derivatives for the advection
		// Set the grid position
		gridPosition[0] = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];
		for (auto l = 0; l < advectionHandlers.size(); l++) {
			advectionHandlers[l]->computePartialsForAdvection(network,
				advecVals.data(), advecIndices.data(), gridPosition, hxLeft,
				hxRight, xi - localXS);

			auto hAdvecVals = HostUnmanaged(advecVals.data(), 2 * nAdvec);
			deep_copy(
				subview(vals, std::make_pair(valIndex, valIndex + 2 * nAdvec)),
				hAdvecVals);
			valIndex += 2 * nAdvec;
		}

		// Get the concentrations at this grid point
		auto concOffset = subview(concs, xi, Kokkos::ALL).view();

		auto surfacePos = grid[surfacePosition + 1];
		auto curXPos = (grid[xi] + grid[xi + 1]) / 2.0;
		auto prevXPos = (grid[xi - 1] + grid[xi]) / 2.0;
		auto curDepth = curXPos - surfacePos;
		auto curSpacing = curXPos - prevXPos;

		// Compute all the partial derivatives for the reactions
		partialDerivativeCounter->increment();
		partialDerivativeTimer->start();
		network.computeAllPartials(concOffset,
			subview(vals, std::make_pair(valIndex, valIndex + nNetworkEntries)),
			xi + 1 - localXS, curDepth, curSpacing);
		partialDerivativeTimer->stop();
		valIndex += nNetworkEntries;
	}
	Kokkos::fence();
	ierr = MatSetValuesCOO(J, vals.data(), ADD_VALUES);
	checkPetscError(
		ierr, "PetscSolver1DHandler::computeJacobian: MatSetValuesCOO failed.");

	// Reset the values
    resetJacobianValues();

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreKokkosOffsetViewDOF(da, localC, &concs);
	checkPetscError(ierr,
		"PetscSolver1DHandler::computeJacobian: "
		"DMDAVecRestoreKokkosOffsetViewDOF failed.");
}

} /* end namespace handler */
} // namespace solver
} // namespace xolotl
