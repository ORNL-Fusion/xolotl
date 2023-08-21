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
	// Degrees of freedom is the total number of clusters in the network
	// + moments
	const auto dof = network.getDOF();

	// We can update the surface position
	// if we are using a restart file
	if (not networkName.empty() and surfaceOffset == 0) {
		io::XFile xfile(networkName);
		auto concGroup = xfile.getGroup<io::XFile::ConcentrationGroup>();
		if (concGroup and concGroup->hasTimesteps()) {
			auto tsGroup = concGroup->getLastTimestepGroup();
			assert(tsGroup);
			grid = tsGroup->readGrid();
		}
	}
	else {
		// Generate the grid in the x direction which will give us the size of
		// the DMDA
		generateGrid(surfaceOffset);
	}

	// Update the number of grid points from the previous loop
	nX = grid.size() - 2;

	// Generate the separate grid for the temperature
	generateTemperatureGrid();

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
		if (isRobin)
			ss << " and Robin for temperature";
		for (auto pair : initialConc) {
			ss << ", initial concentration for Id: " << pair.first
			   << " of: " << pair.second << " nm-3";
		}
		ss << ", grid (nm): ";
		for (auto i = 1; i < grid.size() - 1; i++) {
			ss << grid[i] - grid[1] << " ";
		}
		ss << std::endl;

		if (not sameTemperatureGrid) {
			ss << "Temperature grid (nm): ";
			for (auto i = 0; i < temperatureGrid.size(); i++) {
				ss << temperatureGrid[i] - temperatureGrid[1] << " ";
			}
			ss << std::endl;
		}
		XOLOTL_LOG << ss.str();
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	if (isMirror) {
		PetscCallVoid(DMDACreate1d(
			xolotlComm, DM_BOUNDARY_MIRROR, nX, dof + 1, 1, NULL, &da));
	}
	else {
		PetscCallVoid(DMDACreate1d(
			xolotlComm, DM_BOUNDARY_PERIODIC, nX, dof + 1, 1, NULL, &da));
	}
	PetscCallVoid(DMSetFromOptions(da));
	PetscCallVoid(DMSetUp(da));
}

void
PetscSolver1DHandler::initializeSolverContext(DM& da, TS& ts)
{
	// Degrees of freedom is the total number of clusters in the network
	// + moments
	const auto dof = network.getDOF();

	// Initialize the surface of the first advection handler corresponding to
	// the advection toward the surface (or a dummy one if it is deactivated)
	advectionHandlers[0]->setLocation(0.0);

	/* The ofill (thought of as a dof by dof 2d (row-oriented) array represents
	 * the nonzero coupling between degrees of freedom at one point with
	 * degrees of freedom on the adjacent point to the left or right.
	 */
	std::vector<core::RowColPair> soretEntries;
	std::vector<core::RowColPair> difEntries;
	std::vector<std::vector<core::RowColPair>> advEntries(
		advectionHandlers.size(), std::vector<core::RowColPair>{});

	// Initialize the temperature handler
	temperatureHandler->initialize(dof);

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
	PetscCallVoid(DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL));
	// Set it in the handler
	setLocalCoordinates(xs, xm);

	// Tell the network the number of grid points on this process with ghosts
	// TODO: do we need the ghost points?
	network.setGridSize(localXM + 2);

	// Get the diagonal fill
	network.getDiagonalFill(dfill);

	// The soret initialization needs to be done after the network
	// because it adds connectivities the network would remove
	soretDiffusionHandler->initialize(network, soretEntries, grid, localXS);

	// Load up the block fills
	Mat J;
	PetscCallVoid(TSGetRHSJacobian(ts, &J, nullptr, nullptr, nullptr));
	auto nwEntries = convertToRowColPairList(dof, dfill);
	nNetworkEntries = nwEntries.size();
	//
	// "+ 1" for temperature
	auto dSize = localXM *
		(nNetworkEntries + soretEntries.size() + difEntries.size() + 1);
	// FIXME
	int nAdvec = 0;
	for (auto&& handler : advectionHandlers) {
		nAdvec = std::max(nAdvec, handler->getNumberOfAdvecting());
	}
	auto oSize = localXM * 2 *
		(soretEntries.size() + difEntries.size() + advEntries.size() * nAdvec +
			1);
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
		// Soret diffusion
		for (auto&& component : soretEntries) {
			mapMatStencilsToCoords(component, i, i, rows, cols);
			mapMatStencilsToCoords(component, i, i - 1, rows, cols);
			mapMatStencilsToCoords(component, i, i + 1, rows, cols);
		}
		partialsCount += soretEntries.size() * 3;
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

	nPartials = rows.size();
	PetscCallVoid(
		MatSetPreallocationCOO(J, nPartials, rows.data(), cols.data()));

	// Initialize the arrays for the reaction partial derivatives
	vals = Kokkos::View<double*>("solverPartials", nPartials + 1);

	// Set the size of the partial derivatives vectors
	reactingPartialsForCluster.resize(dof, 0.0);

	// Initialize the flux handler
	fluxHandler->initializeFluxHandler(network, 0, grid);
}

void
PetscSolver1DHandler::initializeConcentration(
	DM& da, Vec& C, DM& oldDA, Vec& oldC)
{
	temperature.clear();

	// Initialize the last temperature at each grid point on this process
	for (auto i = 0; i < localXM + 2; i++) {
		temperature.push_back(0.0);
	}

	// Initialize the grid for the diffusion
	diffusionHandler->initializeDiffusionGrid(
		advectionHandlers, grid, localXM, localXS);
	soretDiffusionHandler->updateSurfacePosition(0);
	temperatureHandler->updateSurfacePosition(0, temperatureGrid);

	// Initialize the grid for the advection
	advectionHandlers[0]->initializeAdvectionGrid(
		advectionHandlers, grid, localXM, localXS);

	// Degrees of freedom is the total number of clusters in the network
	// + moments
	const auto dof = network.getDOF();

	// If this is the first solver loop
	if (surfaceOffset == 0) {
		// Pointer for the concentration vector
		PetscScalar** concentrations = nullptr;
		PetscCallVoid(DMDAVecGetArrayDOF(da, C, &concentrations));

		// Get the last time step written in the HDF5 file
		bool hasConcentrations = false;
		std::unique_ptr<io::XFile> xfile;
		std::unique_ptr<io::XFile::ConcentrationGroup> concGroup;
		if (not networkName.empty()) {
			xfile = std::make_unique<io::XFile>(networkName);
			concGroup = xfile->getGroup<io::XFile::ConcentrationGroup>();
			hasConcentrations = (concGroup and concGroup->hasTimesteps());
		}

		// Pointer for the concentration vector at a specific grid point
		PetscScalar* concOffset = nullptr;

		// Loop on all the grid points
		for (auto i = (PetscInt)localXS - 1;
			 i <= (PetscInt)localXS + (PetscInt)localXM; i++) {
			// Temperature
			plsm::SpaceVector<double, 3> gridPosition{0.0, 0.0, 0.0};
			if (i < 0)
				gridPosition[0] = (temperatureGrid[0] - temperatureGrid[1]) /
					(temperatureGrid[temperatureGrid.size() - 1] -
						temperatureGrid[1]);
			else
				gridPosition[0] =
					((temperatureGrid[i] + temperatureGrid[i + 1]) / 2.0 -
						temperatureGrid[1]) /
					(temperatureGrid[temperatureGrid.size() - 1] -
						temperatureGrid[1]);
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
			if (i >= leftOffset and not hasConcentrations and
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
			auto myConcs =
				tsGroup->readConcentrations(*xfile, localXS, localXM);

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
		auto networkTemp = interpolateTemperature();
		std::vector<double> depths;
		for (auto i = 0; i < networkTemp.size(); i++) {
			if (localXS + i == nX + 1)
				depths.push_back(grid[localXS + i] - grid[1]);
			else
				depths.push_back(
					(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
					grid[1]);
		}
		network.setTemperatures(networkTemp, depths);

		/*
		 Restore vectors
		 */
		PetscCallVoid(DMDAVecRestoreArrayDOF(da, C, &concentrations));
	}
	// Read from the previous vector
	else {
		// Get the boundaries of the old DMDA
		PetscInt oldXs, oldXm;
		PetscCallVoid(
			DMDAGetCorners(oldDA, &oldXs, NULL, NULL, &oldXm, NULL, NULL));

		// Pointers to the PETSc arrays that start at the beginning (xs) of the
		// local array
		PetscScalar **concs = nullptr, **oldConcs = nullptr;
		// Get pointers to vector data
		PetscCallVoid(DMDAVecGetArrayDOFRead(da, C, &concs));
		PetscCallVoid(DMDAVecGetArrayDOF(oldDA, oldC, &oldConcs));

		// Get the procId
		int procId;
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);

		// We have to interpolate between grid points because the grid spacing
		// is changing Loop on the current grid
		for (int xi = 1; xi < nX; xi++) {
			// Compute its distance from the bottom
			double distance = grid[grid.size() - 2] - grid[xi + 1];
			// Loop on the old grid to find the same distance
			for (int i = 1; i < oldGrid.size() - 1; i++) {
				double left = oldGrid[oldGrid.size() - 2] - oldGrid[i];
				double right = oldGrid[oldGrid.size() - 2] - oldGrid[i + 1];
				// Check the distance
				if (distance > right - 1.0e-4) {
					// Create the arrays to receive the data
					PetscScalar *rightConc, *leftConc;

					// Check where all the needed data is located
					int procs[3] = {0, 0, 0};
					if (i - 1 >= oldXs && i - 1 < oldXs + oldXm) {
						procs[0] = procId;
					}
					if (i >= oldXs && i < oldXs + oldXm) {
						procs[1] = procId;
					}
					// Take care of the receive proc
					if (xi >= localXS && xi < localXS + localXM) {
						procs[2] = procId;
					}
					// Get which processor will send and receive the information
					int totalProcs[3] = {0, 0, 0};
					MPI_Allreduce(&procs, &totalProcs, 3, MPI_INT, MPI_SUM,
						MPI_COMM_WORLD);

					// If the left data shares the same process as the new one
					if (totalProcs[0] == totalProcs[2]) {
						if (procId == totalProcs[2]) {
							leftConc = oldConcs[i - 1];
						}
					}
					else {
						// We have to send the data
						// Send the left data
						if (procId == totalProcs[0]) {
							// Send the values
							MPI_Send(&oldConcs[i - 1][0], dof + 1, MPI_DOUBLE,
								totalProcs[2], 2, MPI_COMM_WORLD);
						}
						// Receive the data on the new proc
						if (procId == totalProcs[2]) {
							// Receive the data
							leftConc = new PetscScalar[dof + 1];
							MPI_Recv(leftConc, dof + 1, MPI_DOUBLE,
								totalProcs[0], 2, MPI_COMM_WORLD,
								MPI_STATUS_IGNORE);
						}
					}

					// If the right data shares the same process as the new one
					if (totalProcs[1] == totalProcs[2]) {
						if (procId == totalProcs[2]) {
							rightConc = oldConcs[i];
						}
					}
					else {
						// We have to send the data
						// Send the right data
						if (procId == totalProcs[1]) {
							// Send the values
							MPI_Send(&oldConcs[i][0], dof + 1, MPI_DOUBLE,
								totalProcs[2], 1, MPI_COMM_WORLD);
						}
						// Receive the data on the new proc
						if (procId == totalProcs[2]) {
							// Receive the data
							rightConc = new PetscScalar[dof + 1];
							MPI_Recv(rightConc, dof + 1, MPI_DOUBLE,
								totalProcs[1], 1, MPI_COMM_WORLD,
								MPI_STATUS_IGNORE);
						}
					}

					// Compute the new value on the new proc
					if (procId == totalProcs[2]) {
						// Compute the location of the new grid point within the
						// old segment
						double xFactor = (distance - left) / (right - left);
						// Get the pointer to the data we want to update
						PetscScalar* newConc = concs[xi];
						// Loop on the DOF
						for (int k = 0; k < dof + 1; k++) {
							newConc[k] = leftConc[k] +
								(rightConc[k] - leftConc[k]) * xFactor;
						}

						if (totalProcs[2] != totalProcs[0])
							delete leftConc;
						if (totalProcs[2] != totalProcs[1])
							delete rightConc;
					}

					break;
				}
			}
		}

		// Update the temperature
		// Pointer for the concentration vector at a specific grid point
		PetscScalar* concOffset = nullptr;
		for (auto i = (PetscInt)localXS;
			 i < (PetscInt)localXS + (PetscInt)localXM; i++) {
			concOffset = concs[i];
			temperature[i - localXS + 1] = concOffset[dof];
		}
		if (surfaceOffset > 0 and localXS == 0) {
			temperature[1] = temperature[2];
			concs[0][dof] = temperature[1];
			for (auto pair : initialConc) {
				concs[1][pair.first] = pair.second;
			}
		}
		temperature[0] = temperature[1];
		temperature[localXM + 1] = temperature[localXM];

		// Update the network with the temperature
		auto networkTemp = interpolateTemperature();
		std::vector<double> depths;
		for (auto i = 0; i < networkTemp.size(); i++) {
			if (localXS + i == nX + 1)
				depths.push_back(grid[localXS + i] - grid[1]);
			else
				depths.push_back(
					(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
					grid[1]);
		}
		network.setTemperatures(networkTemp, depths);

		// Restore the vectors
		PetscCallVoid(DMDAVecRestoreArrayDOFRead(da, C, &concs));
		PetscCallVoid(DMDAVecRestoreArrayDOF(oldDA, oldC, &oldConcs));

		// Boundary conditions
		// Set the index to scatter at the surface
		PetscInt *lidxFrom, *lidxTo, lict = 0;
		PetscCallVoid(PetscMalloc1(1, &lidxTo));
		PetscCallVoid(PetscMalloc1(1, &lidxFrom));
		lidxTo[0] = 0;
		lidxFrom[0] = 0;

		// Create the scatter object
		VecScatter scatter;
		IS isTo, isFrom;
		PetscCallVoid(ISCreateBlock(PetscObjectComm((PetscObject)da), dof + 1,
			1, lidxTo, PETSC_OWN_POINTER, &isTo));
		PetscCallVoid(ISCreateBlock(PetscObjectComm((PetscObject)oldDA),
			dof + 1, 1, lidxFrom, PETSC_OWN_POINTER, &isFrom));

		// Create the scatter object
		PetscCallVoid(VecScatterCreate(oldC, isFrom, C, isTo, &scatter));

		// Do the scatter
		PetscCallVoid(
			VecScatterBegin(scatter, oldC, C, INSERT_VALUES, SCATTER_FORWARD));
		PetscCallVoid(
			VecScatterEnd(scatter, oldC, C, INSERT_VALUES, SCATTER_FORWARD));

		// Destroy everything we don't need anymore
		PetscCallVoid(VecScatterDestroy(&scatter));
		PetscCallVoid(ISDestroy(&isTo));
		PetscCallVoid(ISDestroy(&isFrom));

		// Reset the offset
		surfaceOffset = 0;

		//			VecView(oldC, PETSC_VIEWER_STDOUT_WORLD);
		//			VecView(C, PETSC_VIEWER_STDOUT_WORLD);

		// Destroy everything we don't need anymore
		PetscCallVoid(VecDestroy(&oldC));
		PetscCallVoid(DMDestroy(&oldDA));
	}

	return;
}

void
PetscSolver1DHandler::initGBLocation(DM& da, Vec& C)
{
	// Pointer for the concentration vector
	PetscScalar** concentrations = nullptr;
	PetscCallVoid(DMDAVecGetArrayDOF(da, C, &concentrations));

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
	PetscCallVoid(DMDAVecRestoreArrayDOF(da, C, &concentrations));

	return;
}

std::vector<std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
PetscSolver1DHandler::getConcVector(DM& da, Vec& C)
{
	// Initial declaration
	const double* gridPointSolution = nullptr;

	// Pointer for the concentration vector
	PetscScalar** concentrations = nullptr;
	PetscCallContinue(DMDAVecGetArrayDOFRead(da, C, &concentrations));

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
	PetscCallContinue(DMDAVecRestoreArrayDOFRead(da, C, &concentrations));

	return toReturn;
}

void
PetscSolver1DHandler::setConcVector(DM& da, Vec& C,
	std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>&
		concVector)
{
	// Pointer for the concentration vector
	PetscScalar* gridPointSolution = nullptr;
	PetscScalar** concentrations = nullptr;
	PetscCallVoid(DMDAVecGetArrayDOF(da, C, &concentrations));

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
	PetscCallVoid(DMDAVecRestoreArrayDOF(da, C, &concentrations));

	// Get the complete data array, including ghost cells to set the temperature
	// at the ghost points
	Vec localSolution;
	PetscCallVoid(DMGetLocalVector(da, &localSolution));
	PetscCallVoid(DMGlobalToLocalBegin(da, C, INSERT_VALUES, localSolution));
	PetscCallVoid(DMGlobalToLocalEnd(da, C, INSERT_VALUES, localSolution));
	// Get the array of concentration
	PetscCallVoid(DMDAVecGetArrayDOFRead(da, localSolution, &concentrations));

	// Getthe DOF of the network
	const auto dof = network.getDOF();

	// Loop on the grid points
	for (auto i = -1; i <= localXM; ++i) {
		gridPointSolution = concentrations[localXS + i];

		// Get the temperature
		temperature[i + 1] = gridPointSolution[dof];
	}
	// Update the network with the temperature
	auto networkTemp = interpolateTemperature();
	std::vector<double> depths;
	for (auto i = 0; i < networkTemp.size(); i++) {
		if (localXS + i == nX + 1)
			depths.push_back(grid[localXS + i] - grid[1]);
		else
			depths.push_back(
				(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 - grid[1]);
	}
	network.setTemperatures(networkTemp, depths);

	// Restore the solutionArray
	PetscCallVoid(
		DMDAVecRestoreArrayDOFRead(da, localSolution, &concentrations));
	PetscCallVoid(DMRestoreLocalVector(da, &localSolution));

	return;
}

void
PetscSolver1DHandler::updateConcentration(
	TS& ts, Vec& localC, Vec& F, PetscReal ftime)
{
	// Get the local data vector from PETSc
	DM da;
	PetscCallVoid(TSGetDM(ts, &da));

	// Pointers to the PETSc arrays that start at the beginning (localXS) of the
	// local array!
	PetscOffsetView<const PetscScalar**> concs;
	PetscCallVoid(DMDAVecGetKokkosOffsetViewDOF(da, localC, &concs));
	PetscOffsetView<PetscScalar**> updatedConcs;
	PetscCallVoid(DMDAVecGetKokkosOffsetViewDOFWrite(da, F, &updatedConcs));

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
			if (xi < leftOffset || xi > nX - 1 - rightOffset)
				continue;

			// We are only interested in the helium near the surface
			if ((grid[xi] + grid[xi + 1]) / 2.0 - grid[1] > 2.0)
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
		if ((xi == 0 || (xi == nX - 1 && isRobin)) && xi >= localXS &&
			xi < localXS + localXM) {
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

			temperatureHandler->computeTemperature(ftime, concVector.data(),
				updatedConcOffset, hxLeft, hxRight, xi);
		}

		// Compute the old and new array offsets
		auto concOffset = subview(concs, xi, Kokkos::ALL).view();
		auto updatedConcOffset = subview(updatedConcs, xi, Kokkos::ALL).view();

		// Set the grid fraction
		if (xi < 0) {
			gridPosition[0] = (grid[0] - grid[1]) / (grid.back() - grid[1]);
		}
		else {
			gridPosition[0] = ((grid[xi] + grid[xi + 1]) / 2.0 - grid[1]) /
				(grid.back() - grid[1]);
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
		if (xi < leftOffset || xi > nX - 1 - rightOffset) {
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

		// Compute the temperature over the locally owned part of the grid
		if (xi >= localXS && xi < localXS + localXM) {
			temperatureHandler->computeTemperature(ftime, concVector.data(),
				updatedConcOffset, hxLeft, hxRight, xi);
		}
	}

	// Share the information with all the processes
	bool totalTempHasChanged = false;
	auto xolotlComm = util::getMPIComm();
	MPI_Allreduce(&tempHasChanged, &totalTempHasChanged, 1, MPI_C_BOOL, MPI_LOR,
		xolotlComm);

	if (totalTempHasChanged) {
		// Update the network with the temperature
		auto networkTemp = interpolateTemperature();
		std::vector<double> depths;
		for (auto i = 0; i < networkTemp.size(); i++) {
			if (localXS + i == nX + 1)
				depths.push_back(grid[localXS + i] - grid[1]);
			else
				depths.push_back(
					(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
					grid[1]);
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

		// Everything to the left of the surface is empty
		if (xi < leftOffset || xi > nX - 1 - rightOffset) {
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

		if (xi == 0)
			continue;

		// ---- Compute Soret diffusion over the locally owned part of the grid
		// -----
		soretDiffusionHandler->computeDiffusion(network,
			core::StencilConcArray{concVector.data(), 3}, updatedConcOffset,
			hxLeft, hxRight, xi - localXS);

		// ----- Account for flux of incoming particles -----
		fluxHandler->computeIncidentFlux(ftime, updatedConcOffset, xi, 0);

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

		auto surfacePos = grid[1];
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
	PetscCallVoid(DMDAVecRestoreKokkosOffsetViewDOF(da, localC, &concs));
	PetscCallVoid(DMDAVecRestoreKokkosOffsetViewDOFWrite(da, F, &updatedConcs));
}

void
PetscSolver1DHandler::computeJacobian(
	TS& ts, Vec& localC, Mat& J, PetscReal ftime)
{
	// Get the distributed array
	DM da;
	PetscCallVoid(TSGetDM(ts, &da));

	PetscOffsetView<const PetscScalar**> concs;
	PetscCallVoid(DMDAVecGetKokkosOffsetViewDOF(da, localC, &concs));

	// Degrees of freedom is the total number of clusters in the network
	const auto dof = network.getDOF();

	// Get the total number of diffusing clusters
	const auto nSoret =
		std::max(soretDiffusionHandler->getNumberOfDiffusing(), 0);
	const auto nDiff = std::max(diffusionHandler->getNumberOfDiffusing(), 0);

	// Get the total number of advecting clusters
	int nAdvec = 0;
	for (auto&& handler : advectionHandlers) {
		nAdvec = std::max(nAdvec, handler->getNumberOfAdvecting());
	}

	// Arguments for MatSetValuesStencil called below
	PetscScalar tempVals[3];
	IdType tempIndices[1];
	auto diffVals = std::vector<PetscScalar>(3 * nDiff);
	auto diffIndices = std::vector<IdType>(nDiff);
	auto advecVals = std::vector<PetscScalar>(2 * nAdvec);
	auto advecIndices = std::vector<IdType>(nAdvec);
	using ConcSubView = Kokkos::View<const double*>;
	Kokkos::Array<ConcSubView, 3> concVector;
	Kokkos::Array<ConcSubView::host_mirror_type, 3> hConcVec;
	const double* hConcPtrVec[3];
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

		// Get the concentrations at this grid point
		auto concOffset = subview(concs, xi, Kokkos::ALL).view();

		auto tempIndex = valIndex;
		if (xi >= localXS && xi < localXS + localXM) {
			// Fill the concVector with the pointer to the middle, left, and
			// right grid points
			int id = 0;
			for (auto&& xId : {xi, (PetscInt)xi - 1, xi + 1}) {
				concVector[id] = subview(concs, xId, Kokkos::ALL).view();
				hConcVec[id] = create_mirror_view(concVector[id]);
				deep_copy(hConcVec[id], concVector[id]);
				hConcPtrVec[id] = hConcVec[id].data();
				++id;
			}
			valIndex += 3;
		}

		// Heat condition
		if ((xi == 0 || (xi == nX - 1 && isRobin)) && xi >= localXS &&
			xi < localXS + localXM) {
			// Get the partial derivatives for the temperature
			auto setValues = temperatureHandler->computePartialsForTemperature(
				ftime, hConcPtrVec, tempVals, tempIndices, hxLeft, hxRight, xi);

			if (setValues) {
				hTempVals(tempIndex + 0) += tempVals[0];
				hTempVals(tempIndex + 1) += tempVals[1];
				hTempVals(tempIndex + 2) += tempVals[2];
			}
		}

		// Set the grid fraction
		if (xi < 0) {
			gridPosition[0] = (temperatureGrid[0] - temperatureGrid[1]) /
				(temperatureGrid.back() - temperatureGrid[1]);
		}
		else {
			gridPosition[0] =
				((temperatureGrid[xi] + temperatureGrid[xi + 1]) / 2.0 -
					temperatureGrid[1]) /
				(temperatureGrid.back() - temperatureGrid[1]);
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
		if (xi < leftOffset || xi > nX - 1 - rightOffset)
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
				ftime, hConcPtrVec, tempVals, tempIndices, hxLeft, hxRight, xi);

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
		auto networkTemp = interpolateTemperature();
		std::vector<double> depths;
		for (auto i = 0; i < networkTemp.size(); i++) {
			if (localXS + i == nX + 1)
				depths.push_back(grid[localXS + i] - grid[1]);
			else
				depths.push_back(
					(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
					grid[1]);
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
			if (xi < leftOffset || xi > nX - 1 - rightOffset)
				continue;

			// We are only interested in the helium near the surface
			if ((grid[xi] + grid[xi + 1]) / 2.0 - grid[1] > 2.0)
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
		if (xi < leftOffset || xi > nX - 1 - rightOffset) {
			valIndex += 3 * nSoret;
			valIndex += 3 * nDiff;
			valIndex += 2 * nAdvec * advectionHandlers.size();
			valIndex += nNetworkEntries;
			continue;
		}

		// Free surface GB
		if (std::find_if(begin(gbVector), end(gbVector),
				[=](auto&& pair) { return xi == pair[0]; }) != end(gbVector)) {
			// TODO: If the gbVector is initialized before the preallocation, we
			// could simply avoid the extra entries
			valIndex += 3 * nSoret;
			valIndex += 3 * nDiff;
			valIndex += 2 * nAdvec * advectionHandlers.size();
			valIndex += nNetworkEntries;
			continue;
		}

		if (xi == 0)
			continue;

		// Fill the concVector with the pointer to the middle, left, and right
		// grid points
		concVector[0] = subview(concs, xi, Kokkos::ALL).view(); // middle
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

		// Get the partial derivatives for the Soret diffusion
		soretDiffusionHandler->computePartialsForDiffusion(network,
			core::StencilConcArray{concVector.data(), 3},
			subview(vals, std::make_pair(valIndex, valIndex + 3 * nSoret)),
			hxLeft, hxRight, xi - localXS);
		valIndex += 3 * nSoret;

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

		auto surfacePos = grid[1];
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
	PetscCallVoid(MatSetValuesCOO(J, vals.data(), ADD_VALUES));

	// Reset the values
	resetJacobianValues();

	/*
	 Restore vectors
	 */
	PetscCallVoid(DMDAVecRestoreKokkosOffsetViewDOF(da, localC, &concs));
}

} /* end namespace handler */
} // namespace solver
} // namespace xolotl
