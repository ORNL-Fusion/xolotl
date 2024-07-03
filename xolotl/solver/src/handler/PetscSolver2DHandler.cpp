#include <petscconf.h>
#include <petscdmda_kokkos.hpp>

#include <xolotl/core/Constants.h>
#include <xolotl/core/Types.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/io/XFile.h>
#include <xolotl/solver/handler/PetscSolver2DHandler.h>
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
PetscSolver2DHandler::createSolverContext(DM& da)
{
	// Degrees of freedom is the total number of clusters in the network
	// + moments
	const auto dof = network.getDOF();

	bool restarting = this->checkForRestart();

	// We can update the surface position
	// if we are using a restart file
	if (restarting and surfaceOffset == 0) {
		io::XFile xfile(restartFile);
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

	// Set the position of the surface
	for (auto j = 0; j < nY; j++) {
		surfacePosition.push_back(0);
	}

	// We can update the surface position
	// if we are using a restart file
	if (restarting and movingSurface) {
		io::XFile xfile(restartFile);
		auto concGroup = xfile.getGroup<io::XFile::ConcentrationGroup>();
		if (concGroup and concGroup->hasTimesteps()) {
			auto tsGroup = concGroup->getLastTimestepGroup();
			assert(tsGroup);
			auto surfaceIndices = tsGroup->readSurface2D();

			// Set the actual surface positions
			for (auto i = 0; i < surfaceIndices.size(); i++) {
				surfacePosition[i] = surfaceIndices[i];
			}
		}
	}

	// Prints info on one process
	auto xolotlComm = util::getMPIComm();
	int procId = util::getMPIRank();
	if (procId == 0) {
		util::StringStream ss;
		ss << "SolverHandler: 2D simulation with surface BC: ";
		std::string bcString = "periodic";
		if (isMirror)
			bcString = "mirror";
		if (leftOffset == 1)
			ss << "free surface";
		else
			ss << bcString;
		ss << ", bulk BC: ";
		if (rightOffset == 1)
			ss << "free surface";
		else
			ss << bcString;
		ss << ", left BC: ";
		if (topOffset == 1)
			ss << "free surface";
		else
			ss << "periodic";
		ss << ", right BC: ";
		if (bottomOffset == 1)
			ss << "free surface";
		else
			ss << "periodic";
		for (auto pair : initialConc) {
			ss << ", initial concentration for Id: " << pair.first
			   << " of: " << pair.second << " nm-3";
		}
		ss << ", grid (nm): ";
		for (auto i = 1; i < grid.size() - 1; i++) {
			ss << grid[i] - grid[surfacePosition[0] + 1] << " ";
		}
		XOLOTL_LOG << ss.str();
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if (isMirror) {
		PetscCallVoid(DMDACreate2d(xolotlComm, DM_BOUNDARY_MIRROR,
			DM_BOUNDARY_PERIODIC, DMDA_STENCIL_STAR, nX, nY, PETSC_DECIDE,
			PETSC_DECIDE, dof + 1, 1, NULL, NULL, &da));
	}
	else {
		PetscCallVoid(DMDACreate2d(xolotlComm, DM_BOUNDARY_PERIODIC,
			DM_BOUNDARY_PERIODIC, DMDA_STENCIL_STAR, nX, nY, PETSC_DECIDE,
			PETSC_DECIDE, dof + 1, 1, NULL, NULL, &da));
	}
	PetscCallVoid(DMSetFromOptions(da));
	PetscCallVoid(DMSetUp(da));
}

void
PetscSolver2DHandler::initializeSolverContext(DM& da, Mat& J)
{
	// + moments
	const auto dof = network.getDOF();

	// Initialize the surface of the first advection handler corresponding to
	// the advection toward the surface
	advectionHandlers[0]->setLocation(grid[surfacePosition[0] + 1] - grid[1]);

	/* The ofill (thought of as a dof by dof 2d (row-oriented) array represents
	 * the nonzero coupling between degrees of freedom at one point with
	 * degrees of freedom on the adjacent point.
	 */
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
	for (auto i = 0; i < advectionHandlers.size(); i++) {
		advectionHandlers[i]->initialize(network, advEntries[i]);
	}

	// Get the local boundaries
	PetscInt xs, xm, ys, ym;
	PetscCallVoid(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));
	// Set it in the handler
	setLocalCoordinates(xs, xm, ys, ym);

	// Tell the network the number of grid points on this process with ghosts
	// TODO: do we need the ghost points?
	network.setGridSize(localXM + 2);

	// Get the diagonal fill
	network.getDiagonalFill(dfill);

	// Load up the block fills
	auto nwEntries = convertToRowColPairList(dof, dfill);
	nNetworkEntries = nwEntries.size();
	//
	// "+ 1" for temperature
	auto dSize = localYM * localXM * (nNetworkEntries + difEntries.size() + 1);
	int nAdvec = 0;
	for (auto&& handler : advectionHandlers) {
		nAdvec = std::max(nAdvec, handler->getNumberOfAdvecting());
	}
	auto oSize = localYM * localXM * 2 *
		(difEntries.size() + advEntries.size() * nAdvec + 1);
	auto nPartials = dSize + oSize;
	// TODO: should the count used for this reservation be more exact?
	std::vector<PetscInt> rows, cols;
	rows.reserve(nPartials);
	cols.reserve(nPartials);
	auto mapMatStencilsToCoords =
		[da](const core::RowColPair& component, PetscInt yRow, PetscInt xRow,
			PetscInt yCol, PetscInt xCol, std::vector<PetscInt>& gRows,
			std::vector<PetscInt>& gCols) {
			PetscFunctionBeginUser;
			MatStencil stCrds[2];
			stCrds[0].j = yRow;
			stCrds[0].i = xRow;
			stCrds[0].c = component[0];
			stCrds[1].j = yCol;
			stCrds[1].i = xCol;
			stCrds[1].c = component[1];
			PetscInt coo[2];
			PetscCall(DMDAMapMatStencilToGlobal(da, 2, stCrds, coo));
			gRows.push_back(coo[0]);
			gCols.push_back(coo[1]);
			PetscFunctionReturn(0);
		};
	std::size_t partialsCount{0};
	for (auto j = localYS; j < localYS + localYM; ++j) {
		for (auto i = localXS; i < localXS + localXM; ++i) {
			// temperature
			mapMatStencilsToCoords({dof, dof}, j, i, j, i, rows, cols);
			mapMatStencilsToCoords({dof, dof}, j, i, j, i - 1, rows, cols);
			mapMatStencilsToCoords({dof, dof}, j, i, j, i + 1, rows, cols);
			mapMatStencilsToCoords({dof, dof}, j, i, j - 1, i, rows, cols);
			mapMatStencilsToCoords({dof, dof}, j, i, j + 1, i, rows, cols);
			partialsCount += 5;
		}
	}
	for (auto j = localYS; j < localYS + localYM; ++j) {
		for (auto i = localXS; i < localXS + localXM; ++i) {
			// diffusion
			for (auto&& component : difEntries) {
				mapMatStencilsToCoords(component, j, i, j, i, rows, cols);
				mapMatStencilsToCoords(component, j, i, j, i - 1, rows, cols);
				mapMatStencilsToCoords(component, j, i, j, i + 1, rows, cols);
				mapMatStencilsToCoords(component, j, i, j - 1, i, rows, cols);
				mapMatStencilsToCoords(component, j, i, j + 1, i, rows, cols);
			}
			partialsCount += difEntries.size() * 5;
			// advection
			plsm::SpaceVector<double, 3> pos{0.0, 0.0, 0.0};
			pos[0] = (grid[i] + grid[i + 1]) / 2.0 - grid[1];
			for (std::size_t l = 0; l < advectionHandlers.size(); ++l) {
				auto offsets =
					advectionHandlers[l]->getStencilForAdvection(pos);
				if (advectionHandlers[l]->isPointOnSink(pos)) {
					for (auto&& component : advEntries[l]) {
						mapMatStencilsToCoords(component, j, i, j - offsets[1],
							i - offsets[0], rows, cols);
						mapMatStencilsToCoords(component, j, i, j + offsets[1],
							i + offsets[0], rows, cols);
					}
				}
				else {
					for (auto&& component : advEntries[l]) {
						mapMatStencilsToCoords(
							component, j, i, j, i, rows, cols);
						mapMatStencilsToCoords(component, j, i, j + offsets[1],
							i + offsets[0], rows, cols);
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
				mapMatStencilsToCoords(component, j, i, j, i, rows, cols);
			}
			partialsCount += nNetworkEntries;
		}
	}
	nPartials = rows.size();
	PetscCallVoid(
		MatSetPreallocationCOO(J, nPartials, rows.data(), cols.data()));

	// Initialize the arrays for the reaction partial derivatives
	vals = Kokkos::View<double*>("solverPartials", nPartials + 1);

	// Set the size of the partial derivatives vectors
	reactingPartialsForCluster.resize(dof, 0.0);

	// Initialize the flux handler
	fluxHandler->initializeFluxHandler(network, surfacePosition[0], grid);
}

void
PetscSolver2DHandler::initializeConcentration(
	DM& da, Vec& C, DM& oldDA, Vec& oldC)
{
	temperature.clear();

	// Initialize the last temperature at each grid point on this process
	for (auto i = 0; i < localXM + 2; i++) {
		temperature.push_back(0.0);
	}

	// Initialize the grid for the diffusion
	diffusionHandler->initializeDiffusionGrid(
		advectionHandlers, grid, localXM, localXS, localYM, hY, localYS);

	// Initialize the grid for the advection
	advectionHandlers[0]->initializeAdvectionGrid(
		advectionHandlers, grid, localXM, localXS, localYM, hY, localYS);

	// Degrees of freedom is the total number of clusters in the network
	// + moments
	const auto dof = network.getDOF();

	// If this is the first solver loop
	if (surfaceOffset == 0) {
		// Pointer for the concentration vector
		PetscScalar*** concentrations = nullptr;
		PetscCallVoid(DMDAVecGetArrayDOF(da, C, &concentrations));

		// Get the last time step written in the HDF5 file
		bool hasConcentrations = false;
		std::unique_ptr<io::XFile> xfile;
		std::unique_ptr<io::XFile::ConcentrationGroup> concGroup;
		if (this->checkForRestart()) {
			xfile = std::make_unique<io::XFile>(restartFile);
			concGroup = xfile->getGroup<io::XFile::ConcentrationGroup>();
			hasConcentrations = (concGroup and concGroup->hasTimesteps());
		}

		// Pointer for the concentration vector at a specific grid point
		PetscScalar* concOffset = nullptr;

		// Loop on all the grid points
		for (auto j = localYS; j < localYS + localYM; j++) {
			for (auto i = (PetscInt)localXS - 1;
				 i <= (PetscInt)localXS + (PetscInt)localXM; i++) {
				// Temperature
				plsm::SpaceVector<double, 3> gridPosition{0.0, j * hY, 0.0};
				if (i < 0)
					gridPosition[0] = (grid[0] - grid[surfacePosition[j] + 1]) /
						(grid[grid.size() - 1] - grid[surfacePosition[j] + 1]);
				else
					gridPosition[0] = ((grid[i] + grid[i + 1]) / 2.0 -
										  grid[surfacePosition[j] + 1]) /
						(grid[grid.size() - 1] - grid[surfacePosition[j] + 1]);
				auto temp =
					temperatureHandler->getTemperature(gridPosition, 0.0);
				temperature[i - localXS + 1] = temp;

				// Boundary conditions
				if (i < localXS || i >= localXS + localXM)
					continue;

				concOffset = concentrations[j][i];
				concOffset[dof] = temp;

				// Loop on all the clusters to initialize at 0.0
				for (auto n = 0; n < dof; n++) {
					concOffset[n] = 0.0;
				}

				// Initialize the option specified concentration
				if (i >= surfacePosition[j] + leftOffset and
					not hasConcentrations and i < nX - rightOffset and
					j >= bottomOffset and j < nY - topOffset) {
					for (auto pair : initialConc) {
						concOffset[pair.first] = pair.second;
					}
				}
			}
		}

		// If the concentration must be set from the HDF5 file
		if (hasConcentrations) {
			assert(concGroup);
			auto tsGroup = concGroup->getLastTimestepGroup();
			assert(tsGroup);

			// Loop on the full grid
			for (auto j = 0; j < nY; j++) {
				for (auto i = 0; i < nX; i++) {
					// Read the concentrations from the HDF5 file
					auto concVector = tsGroup->readGridPoint(i, j);

					// Change the concentration only if we are on the locally
					// owned part of the grid
					if (i >= localXS && i < localXS + localXM && j >= localYS &&
						j < localYS + localYM) {
						concOffset = concentrations[j][i];
						// Loop on the concVector size
						for (auto l = 0; l < concVector.size(); l++) {
							concOffset[(IdType)concVector.at(l).at(0)] =
								concVector.at(l).at(1);
						}
						// Get the temperature
						double temp =
							concVector.at(concVector.size() - 1).at(1);
						temperature[i - localXS + 1] = temp;
					}
				}
			}
		}

		// Update the network with the temperature
		std::vector<double> depths;
		for (auto i = 0; i < temperature.size(); i++) {
			if (localXS + i == nX + 1)
				depths.push_back(
					grid[localXS + i] - grid[surfacePosition[localYS] + 1]);
			else
				depths.push_back(
					(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
					grid[surfacePosition[localYS] + 1]);
		}
		network.setTemperatures(temperature, depths);

		/*
		 Restore vectors
		 */
		PetscCallVoid(DMDAVecRestoreArrayDOF(da, C, &concentrations));
	}
	// Read from the previous vector
	else {
		// Get the boundaries of the old DMDA
		PetscInt oldXs, oldXm, oldYs, oldYm;
		PetscCallVoid(
			DMDAGetCorners(oldDA, &oldXs, &oldYs, NULL, &oldXm, &oldYm, NULL));

		// Pointers to the PETSc arrays that start at the beginning (xs, ys) of
		// the local array
		PetscScalar ***concs = nullptr, ***oldConcs = nullptr;
		// Get pointers to vector data
		PetscCallVoid(DMDAVecGetArrayDOFRead(da, C, &concs));
		PetscCallVoid(DMDAVecGetArrayDOF(oldDA, oldC, &oldConcs));

		// Get the procId
		int procId;
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);

		// We have to interpolate between grid points because the grid spacing
		// is changing
		for (int yj = 1; yj < nY; yj++) {
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
						std::vector<PetscScalar> leftConcVec, rightConcVec;
						PetscScalar *rightConc, *leftConc;

						// Check where all the needed data is located
						int procs[3] = {0, 0, 0};
						if (i - 1 >= oldXs and i - 1 < oldXs + oldXm and
							yj >= oldYs and yj < oldYs + oldYm) {
							procs[0] = procId;
						}
						if (i >= oldXs and i < oldXs + oldXm and yj >= oldYs and
							yj < oldYs + oldYm) {
							procs[1] = procId;
						}
						// Take care of the receive proc
						if (xi >= localXS and xi < localXS + localXM and
							yj >= localYS and yj < localYS + localYM) {
							procs[2] = procId;
						}
						// Get which processor will send and receive the
						// information
						int totalProcs[3] = {0, 0, 0};
						MPI_Allreduce(&procs, &totalProcs, 3, MPI_INT, MPI_SUM,
							MPI_COMM_WORLD);

						// If the left data shares the same process as the new
						// one
						if (totalProcs[0] == totalProcs[2]) {
							if (procId == totalProcs[2]) {
								leftConc = oldConcs[yj][i - 1];
							}
						}
						else {
							// We have to send the data
							// Send the left data
							if (procId == totalProcs[0]) {
								// Send the values
								MPI_Send(&oldConcs[yj][i - 1][0], dof + 1,
									MPI_DOUBLE, totalProcs[2], 2,
									MPI_COMM_WORLD);
							}
							// Receive the data on the new proc
							if (procId == totalProcs[2]) {
								// Receive the data
								leftConcVec.resize(dof + 1);
								leftConc = leftConcVec.data();
								MPI_Recv(leftConc, dof + 1, MPI_DOUBLE,
									totalProcs[0], 2, MPI_COMM_WORLD,
									MPI_STATUS_IGNORE);
							}
						}

						// If the right data shares the same process as the new
						// one
						if (totalProcs[1] == totalProcs[2]) {
							if (procId == totalProcs[2]) {
								rightConc = oldConcs[yj][i];
							}
						}
						else {
							// We have to send the data
							// Send the right data
							if (procId == totalProcs[1]) {
								// Send the values
								MPI_Send(&oldConcs[yj][i][0], dof + 1,
									MPI_DOUBLE, totalProcs[2], 1,
									MPI_COMM_WORLD);
							}
							// Receive the data on the new proc
							if (procId == totalProcs[2]) {
								// Receive the data
								rightConcVec.resize(dof + 1);
								rightConc = rightConcVec.data();
								MPI_Recv(rightConc, dof + 1, MPI_DOUBLE,
									totalProcs[1], 1, MPI_COMM_WORLD,
									MPI_STATUS_IGNORE);
							}
						}

						// Compute the new value on the new proc
						if (procId == totalProcs[2]) {
							// Compute the location of the new grid point within
							// the old segment
							double xFactor = (distance - left) / (right - left);
							// Get the pointer to the data we want to update
							PetscScalar* newConc = concs[yj][xi];
							// Loop on the DOF
							for (int k = 0; k < dof + 1; k++) {
								newConc[k] = leftConc[k] +
									(rightConc[k] - leftConc[k]) * xFactor;
							}
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
				concOffset = concs[yj][i];
				temperature[i - localXS + 1] = concOffset[dof];
			}
			if (surfaceOffset > 0 and localXS == 0) {
				temperature[1] = temperature[2];
				concs[yj][0][dof] = temperature[1];
				for (auto pair : initialConc) {
					concs[yj][1][pair.first] = pair.second;
				}
			}
			temperature[0] = temperature[1];
			temperature[localXM + 1] = temperature[localXM];
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

		// Restore the vectors
		PetscCallVoid(DMDAVecRestoreArrayDOFRead(da, C, &concs));
		PetscCallVoid(DMDAVecRestoreArrayDOF(oldDA, oldC, &oldConcs));

		// Boundary conditions
		// Set the index to scatter at the surface
		PetscInt* lidxTo{nullptr};
		PetscInt* lidxFrom{nullptr};
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

		// Destroy everything we don't need anymore
		PetscCallVoid(VecDestroy(&oldC));
		PetscCallVoid(DMDestroy(&oldDA));
	}

	return;
}

void
PetscSolver2DHandler::initGBLocation(DM& da, Vec& C)
{
	// Pointer for the concentration vector
	PetscScalar*** concentrations = nullptr;
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
		auto yj = std::get<1>(pair);
		// Check if we are on the right process
		if (xi >= localXS && xi < localXS + localXM && yj >= localYS &&
			yj < localYS + localYM) {
			// Get the local concentration
			concOffset = concentrations[yj][xi];

			using HostUnmanaged = Kokkos::View<double*, Kokkos::HostSpace,
				Kokkos::MemoryUnmanaged>;
			auto hConcs = HostUnmanaged(concOffset, dof);
			auto dConcs = Kokkos::View<double*>("Concentrations", dof);
			deep_copy(dConcs, hConcs);

			// Transfer the local amount of Xe clusters
			setLocalXeRate(
				neNetwork.getTotalAtomConcentration(dConcs, Spec::Xe, 1),
				xi - localXS, yj - localYS);

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
PetscSolver2DHandler::getConcVector(DM& da, Vec& C)
{
	// Initial declaration
	const double* gridPointSolution = nullptr;

	// Pointer for the concentration vector
	PetscScalar*** concentrations = nullptr;
	PetscCallContinue(DMDAVecGetArrayDOFRead(da, C, &concentrations));

	// Get the network and dof
	auto& network = getNetwork();
	const auto dof = network.getDOF();

	// Create the vector for the concentrations
	std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
		toReturn;
	std::vector<std::vector<std::vector<std::pair<IdType, double>>>>
		tempTempTempVector;

	// Loop on the grid points
	for (auto j = 0; j < localYM; ++j) {
		std::vector<std::vector<std::pair<IdType, double>>> tempTempVector;
		for (auto i = 0; i < localXM; ++i) {
			gridPointSolution = concentrations[localYS + j][localXS + i];

			// Create the temporary vector for this grid point
			std::vector<std::pair<IdType, double>> tempVector;
			for (auto l = 0; l < dof + 1; ++l) {
				if (std::fabs(gridPointSolution[l]) > 1.0e-16) {
					tempVector.push_back(
						std::make_pair(l, gridPointSolution[l]));
				}
			}
			tempTempVector.push_back(tempVector);
		}
		tempTempTempVector.push_back(tempTempVector);
	}
	toReturn.push_back(tempTempTempVector);

	// Restore the solutionArray
	PetscCallContinue(DMDAVecRestoreArrayDOFRead(da, C, &concentrations));

	return toReturn;
}

void
PetscSolver2DHandler::setConcVector(DM& da, Vec& C,
	std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>&
		concVector)
{
	// Pointer for the concentration vector
	PetscScalar* gridPointSolution = nullptr;
	PetscScalar*** concentrations = nullptr;
	PetscCallVoid(DMDAVecGetArrayDOF(da, C, &concentrations));

	// Loop on the grid points
	for (auto j = 0; j < localYM; ++j) {
		for (auto i = 0; i < localXM; ++i) {
			gridPointSolution = concentrations[localYS + j][localXS + i];

			// Loop on the given vector
			for (auto l = 0; l < concVector[0][j][i].size(); l++) {
				gridPointSolution[concVector[0][j][i][l].first] =
					concVector[0][j][i][l].second;
			}
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
	for (auto j = 0; j < localYM; ++j) {
		for (auto i = -1; i <= localXM; ++i) {
			gridPointSolution = concentrations[localYS + j][localXS + i];

			// Get the temperature
			temperature[i + 1] = gridPointSolution[dof];
		}
		// Update the network with the temperature
		std::vector<double> depths;
		for (auto i = 0; i < temperature.size(); i++) {
			if (localXS + i == nX + 1)
				depths.push_back(
					grid[localXS + i] - grid[surfacePosition[localYS] + 1]);
			else
				depths.push_back(
					(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
					grid[surfacePosition[localYS] + 1]);
		}
		network.setTemperatures(temperature, depths);
	}

	// Restore the solutionArray
	PetscCallVoid(
		DMDAVecRestoreArrayDOFRead(da, localSolution, &concentrations));
	PetscCallVoid(DMRestoreLocalVector(da, &localSolution));

	return;
}

void
PetscSolver2DHandler::updateConcentration(
	TS& ts, Vec& localC, Vec& F, PetscReal ftime)
{
	// Get the local data vector from PETSc
	DM da;
	PetscCallVoid(TSGetDM(ts, &da));

	// Pointers to the PETSc arrays that start at the beginning (localXS,
	// localYS) of the local array
	PetscOffsetView<const PetscScalar***> concs;
	PetscCallVoid(DMDAVecGetKokkosOffsetViewDOF(da, localC, &concs));
	PetscOffsetView<PetscScalar***> updatedConcs;
	PetscCallVoid(DMDAVecGetKokkosOffsetViewDOFWrite(da, F, &updatedConcs));

	// Set some step size variable
	double sy = 1.0 / (hY * hY);

	// Declarations for variables used in the loop
	using ConcSubView = Kokkos::View<const double*>;
	Kokkos::Array<ConcSubView, 5> concVector;
	plsm::SpaceVector<double, 3> gridPosition{0.0, 0.0, 0.0};
	std::vector<double> incidentFluxVector;
	double atomConc = 0.0, totalAtomConc = 0.0;

	// Loop over grid points first for the temperature, including the ghost
	// points in X
	for (auto yj = localYS; yj < localYS + localYM; yj++) {
		temperatureHandler->updateSurfacePosition(surfacePosition[yj], grid);
		bool tempHasChanged = false;
		for (auto xi = (PetscInt)localXS - 1;
			 xi <= (PetscInt)localXS + (PetscInt)localXM; xi++) {
			// Heat condition
			if (xi == surfacePosition[yj] && xi >= localXS &&
				xi < localXS + localXM) {
				// Compute the old and new array offsets
				auto concOffset = subview(concs, yj, xi, Kokkos::ALL).view();
				auto updatedConcOffset =
					subview(updatedConcs, yj, xi, Kokkos::ALL).view();

				// Fill the concVector with the pointer to the middle, left, and
				// right grid points
				concVector[0] = concOffset; // middle
				concVector[1] =
					subview(concs, yj, (PetscInt)xi - 1, Kokkos::ALL)
						.view(); // left
				concVector[2] =
					subview(concs, yj, xi + 1, Kokkos::ALL).view(); // right
				concVector[3] =
					subview(concs, (PetscInt)yj - 1, xi, Kokkos::ALL)
						.view(); // bottom
				concVector[4] =
					subview(concs, yj + 1, xi, Kokkos::ALL).view(); // top

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

				temperatureHandler->computeTemperature(ftime, concVector.data(),
					updatedConcOffset, hxLeft, hxRight, xi, sy, yj);
			}

			// Compute the old and new array offsets
			auto concOffset = subview(concs, yj, xi, Kokkos::ALL).view();
			auto updatedConcOffset =
				subview(updatedConcs, yj, xi, Kokkos::ALL).view();

			// Set the grid fraction
			if (xi < 0)
				gridPosition[0] = (grid[0] - grid[surfacePosition[yj] + 1]) /
					(grid[grid.size() - 1] - grid[surfacePosition[yj] + 1]);
			else
				gridPosition[0] = ((grid[xi] + grid[xi + 1]) / 2.0 -
									  grid[surfacePosition[yj] + 1]) /
					(grid[grid.size() - 1] - grid[surfacePosition[yj] + 1]);
			gridPosition[1] = yj / nY;

			// Get the temperature from the temperature handler
			temperatureHandler->setTemperature(concOffset);
			double temp =
				temperatureHandler->getTemperature(gridPosition, ftime);

			// Update the network if the temperature changed
			if (std::fabs(temperature[xi + 1 - localXS] - temp) > 0.1) {
				temperature[xi + 1 - localXS] = temp;
				tempHasChanged = true;
			}

			// Boundary conditions
			// Everything to the left of the surface is empty
			if (xi < surfacePosition[yj] + leftOffset ||
				xi > nX - 1 - rightOffset) {
				continue;
			}
			// Free surface GB
			bool skip = false;
			for (auto& pair : gbVector) {
				if (xi == std::get<0>(pair) && yj == std::get<1>(pair)) {
					skip = true;
					break;
				}
			}
			if (skip)
				continue;

			// Fill the concVector with the pointer to the middle, left, and
			// right grid points
			concVector[0] = concOffset; // middle
			concVector[1] = subview(concs, yj, (PetscInt)xi - 1, Kokkos::ALL)
								.view(); // left
			concVector[2] =
				subview(concs, yj, xi + 1, Kokkos::ALL).view(); // right
			concVector[3] = subview(concs, (PetscInt)yj - 1, xi, Kokkos::ALL)
								.view(); // bottom
			concVector[4] =
				subview(concs, yj + 1, xi, Kokkos::ALL).view(); // top

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

			// Compute the temperature over the locally owned part of the grid
			if (xi >= localXS && xi < localXS + localXM) {
				temperatureHandler->computeTemperature(ftime, concVector.data(),
					updatedConcOffset, hxLeft, hxRight, xi, sy, yj);
			}
		}

		// TODO: it is updated T more than once per MPI process in preparation
		// of T depending on more than X
		if (tempHasChanged) {
			// Update the network with the temperature
			std::vector<double> depths;
			for (auto i = 0; i < temperature.size(); i++) {
				if (localXS + i == nX + 1)
					depths.push_back(
						grid[localXS + i] - grid[surfacePosition[localYS] + 1]);
				else
					depths.push_back(
						(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
						grid[surfacePosition[localYS] + 1]);
			}
			network.setTemperatures(temperature, depths);
		}
	}

	// Loop over grid points
	for (auto yj = bottomOffset; yj < nY - topOffset; yj++) {
		// Computing the trapped atom concentration is only needed for the
		// attenuation
		if (useAttenuation) {
			// Compute the total concentration of atoms contained in bubbles
			atomConc = 0.0;

			auto& psiNetwork =
				dynamic_cast<core::network::IPSIReactionNetwork&>(network);

			// Loop over grid points
			for (auto xi = surfacePosition[yj] + leftOffset;
				 xi < nX - rightOffset; xi++) {
				// We are only interested in the helium near the surface
				if ((grid[xi] + grid[xi + 1]) / 2.0 -
						grid[surfacePosition[yj] + 1] >
					2.0)
					continue;

				// Check if we are on the right processor
				if (xi >= localXS && xi < localXS + localXM && yj >= localYS &&
					yj < localYS + localYM) {
					// Get the concentrations at this grid point
					auto concOffset =
						subview(concs, yj, xi, Kokkos::ALL).view();

					// Sum the total atom concentration
					atomConc += psiNetwork.getTotalTrappedHeliumConcentration(
									concOffset, 0) *
						(grid[xi + 1] - grid[xi]);
				}
			}

			// Share the concentration with all the processes
			totalAtomConc = 0.0;
			auto xolotlComm = util::getMPIComm();
			MPI_Allreduce(
				&atomConc, &totalAtomConc, 1, MPI_DOUBLE, MPI_SUM, xolotlComm);

			// Set the disappearing rate in the modified TM handler
			psiNetwork.updateTrapMutationDisappearingRate(totalAtomConc);
		}

		// Skip if we are not on the right process
		if (yj < localYS || yj >= localYS + localYM)
			continue;

		// Set the grid position
		gridPosition[1] = yj * hY;

		// Initialize the flux, advection handlers which depend
		// on the surface position at Y
		fluxHandler->initializeFluxHandler(network, surfacePosition[yj], grid);
		advectionHandlers[0]->setLocation(
			grid[surfacePosition[yj] + 1] - grid[1]);

		for (auto xi = localXS; xi < localXS + localXM; xi++) {
			// Compute the old and new array offsets
			auto concOffset = subview(concs, yj, xi, Kokkos::ALL).view();
			auto updatedConcOffset =
				subview(updatedConcs, yj, xi, Kokkos::ALL).view();

			// Fill the concVector with the pointer to the middle, left, right,
			// bottom, and top grid points
			concVector[0] = concOffset; // middle
			concVector[1] = subview(concs, yj, (PetscInt)xi - 1,
				Kokkos::ALL)
								.view(); // left
			concVector[2] =
				subview(concs, yj, xi + 1, Kokkos::ALL).view(); // right
			concVector[3] = subview(concs, (PetscInt)yj - 1, xi,
				Kokkos::ALL)
								.view(); // bottom
			concVector[4] =
				subview(concs, yj + 1, xi, Kokkos::ALL).view(); // top

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
			if (xi < surfacePosition[yj] + leftOffset ||
				xi > nX - 1 - rightOffset || yj < bottomOffset ||
				yj > nY - 1 - topOffset) {
				continue;
			}
			// Free surface GB
			bool skip = false;
			for (auto& pair : gbVector) {
				if (xi == std::get<0>(pair) && yj == std::get<1>(pair)) {
					skip = true;
					break;
				}
			}
			if (skip)
				continue;

			// ----- Account for flux of incoming particles -----
			fluxHandler->computeIncidentFlux(
				ftime, updatedConcOffset, xi, surfacePosition[yj]);

			// ---- Compute diffusion over the locally owned part of the grid
			// -----
			diffusionHandler->computeDiffusion(network,
				core::StencilConcArray{concVector.data(), 5}, updatedConcOffset,
				hxLeft, hxRight, xi - localXS, sy, yj - localYS);

			// ---- Compute advection over the locally owned part of the grid
			// ----- Set the grid position
			gridPosition[0] = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];
			for (auto i = 0; i < advectionHandlers.size(); i++) {
				advectionHandlers[i]->computeAdvection(network, gridPosition,
					core::StencilConcArray{concVector.data(), 5},
					updatedConcOffset, hxLeft, hxRight, xi - localXS, hY,
					yj - localYS);
			}

			auto surfacePos = grid[surfacePosition[yj] + 1];
			auto curXPos = (grid[xi] + grid[xi + 1]) / 2.0;
			auto prevXPos = (grid[xi - 1] + grid[xi]) / 2.0;
			auto curDepth = curXPos - surfacePos;
			auto curSpacing = curXPos - prevXPos;

			// ----- Compute the reaction fluxes over the locally owned part of
			// the grid -----
			fluxCounter->increment();
			fluxTimer->start();
			network.computeAllFluxes(concOffset, updatedConcOffset,
				xi + 1 - localXS, curDepth, curSpacing);
			fluxTimer->stop();
		}
	}

	/*
	 Restore vectors
	 */
	PetscCallVoid(DMDAVecRestoreKokkosOffsetViewDOF(da, localC, &concs));
	PetscCallVoid(DMDAVecRestoreKokkosOffsetViewDOFWrite(da, F, &updatedConcs));
}

void
PetscSolver2DHandler::computeJacobian(
	TS& ts, Vec& localC, Mat& J, PetscReal ftime)
{
	// Get the distributed array
	DM da;
	PetscCallVoid(TSGetDM(ts, &da));

	// Get pointers to vector data
	PetscOffsetView<const PetscScalar***> concs;
	PetscCallVoid(DMDAVecGetKokkosOffsetViewDOF(da, localC, &concs));

	using ConcSubView = Kokkos::View<const double*>;
	Kokkos::Array<ConcSubView, 5> concVector;

	// Get the total number of diffusing clusters
	const auto nDiff = std::max(diffusionHandler->getNumberOfDiffusing(), 0);

	// Get the total number of advecting clusters
	int nAdvec = 0;
	for (auto&& handler : advectionHandlers) {
		nAdvec = std::max(nAdvec, handler->getNumberOfAdvecting());
	}

	// Setup some step size variables
	double sy = 1.0 / (hY * hY);

	// Declarations for variables used in the loop
	double atomConc = 0.0, totalAtomConc = 0.0;
	PetscScalar tempVals[5];
	IdType tempIndices[1];
	Kokkos::Array<ConcSubView::host_mirror_type, 5> hConcVec;
	const double* hConcPtrVec[5];
	plsm::SpaceVector<double, 3> gridPosition{0.0, 0.0, 0.0};

	/*
	 Loop over grid points for the temperature, including ghosts
	 */
	auto hTempVals = Kokkos::View<double*, Kokkos::HostSpace>(
		"Host Temp Jac Vals", localYM * localXM * 5);
	std::size_t valIndex = 0;
	for (PetscInt yj = localYS; yj < localYS + localYM; yj++) {
		temperatureHandler->updateSurfacePosition(surfacePosition[yj], grid);
		bool tempHasChanged = false;
		for (auto xi = (PetscInt)localXS - 1;
			 xi <= (PetscInt)localXS + (PetscInt)localXM; xi++) {
			// Compute the left and right hx
			double hxLeft = 0.0, hxRight = 0.0;
			if (xi >= 1 && xi < nX) {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			}
			else if (xi == -1) {
				hxLeft = 0.0;
				hxRight = (grid[xi + 2] + grid[xi + 1]) / 2.0;
			}
			else if (xi < 1) {
				hxLeft = grid[xi + 1] - grid[xi];
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			}
			else {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = grid[xi + 1] - grid[xi];
			}

			// Get the concentrations at this grid point
			auto concOffset = subview(concs, yj, xi, Kokkos::ALL).view();

			auto tempIndex = valIndex;
			if (xi >= localXS && xi < localXS + localXM) {
				int id = 0;
				using A = std::array<PetscInt, 2>;
				for (auto&& mId : {A{yj, xi}, A{yj, xi - 1}, A{yj, xi + 1},
						 A{yj - 1, xi}, A{yj + 1, xi}}) {
					concVector[id] =
						subview(concs, mId[0], mId[1], Kokkos::ALL).view();
					hConcVec[id] = create_mirror_view(concVector[id]);
					deep_copy(hConcVec[id], concVector[id]);
					hConcPtrVec[id] = hConcVec[id].data();
					++id;
				}
				valIndex += 5;
			}

			// Heat condition
			if (xi == surfacePosition[yj] && xi >= localXS &&
				xi < localXS + localXM) {
				// Get the partial derivatives for the temperature
				auto setValues =
					temperatureHandler->computePartialsForTemperature(ftime,
						hConcPtrVec, tempVals, tempIndices, hxLeft, hxRight, xi,
						sy, yj);

				if (setValues) {
					hTempVals(tempIndex + 0) += tempVals[0];
					hTempVals(tempIndex + 1) += tempVals[1];
					hTempVals(tempIndex + 2) += tempVals[2];
					hTempVals(tempIndex + 3) += tempVals[3];
					hTempVals(tempIndex + 4) += tempVals[4];
				}
			}

			// Set the grid fraction
			if (xi < 0)
				gridPosition[0] = (grid[0] - grid[surfacePosition[yj] + 1]) /
					(grid[grid.size() - 1] - grid[surfacePosition[yj] + 1]);
			else
				gridPosition[0] = ((grid[xi] + grid[xi + 1]) / 2.0 -
									  grid[surfacePosition[yj] + 1]) /
					(grid[grid.size() - 1] - grid[surfacePosition[yj] + 1]);
			gridPosition[1] = yj / nY;

			// Get the temperature from the temperature handler
			temperatureHandler->setTemperature(concOffset);
			double temp =
				temperatureHandler->getTemperature(gridPosition, ftime);

			// Update the network if the temperature changed
			if (std::fabs(temperature[xi + 1 - localXS] - temp) > 0.1) {
				temperature[xi + 1 - localXS] = temp;
				tempHasChanged = true;
			}

			// Boundary conditions
			// Everything to the left of the surface is empty
			if (xi < surfacePosition[yj] + leftOffset ||
				xi > nX - 1 - rightOffset)
				continue;
			// Free surface GB
			if (std::find_if(begin(gbVector), end(gbVector), [=](auto&& pair) {
					return xi == pair[0] && yj == pair[1];
				}) != end(gbVector)) {
				continue;
			}

			// Get the partial derivatives for the temperature
			if (xi >= localXS && xi < localXS + localXM) {
				auto setValues =
					temperatureHandler->computePartialsForTemperature(ftime,
						hConcPtrVec, tempVals, tempIndices, hxLeft, hxRight, xi,
						sy, yj);

				if (setValues) {
					hTempVals(tempIndex + 0) += tempVals[0];
					hTempVals(tempIndex + 1) += tempVals[1];
					hTempVals(tempIndex + 2) += tempVals[2];
					hTempVals(tempIndex + 3) += tempVals[3];
					hTempVals(tempIndex + 4) += tempVals[4];
				}
			}
		}

		if (tempHasChanged) {
			// Update the network with the temperature
			std::vector<double> depths;
			for (auto i = 0; i < temperature.size(); i++) {
				if (localXS + i == nX + 1)
					depths.push_back(
						grid[localXS + i] - grid[surfacePosition[localYS] + 1]);
				else
					depths.push_back(
						(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
						grid[surfacePosition[localYS] + 1]);
			}
			network.setTemperatures(temperature, depths);
		}
	}
	deep_copy(subview(vals, std::make_pair(IdType{0}, localYM * localXM * 5)),
		hTempVals);

	// Loop over the grid points
	for (auto yj = localYS; yj < localYS + localYM; yj++) {
		// Computing the trapped atom concentration is only needed for the
		// attenuation
		if (useAttenuation) {
			// Compute the total concentration of atoms contained in bubbles
			atomConc = 0.0;

			auto& psiNetwork =
				dynamic_cast<core::network::IPSIReactionNetwork&>(network);

			// Loop over grid points
			for (auto xi = surfacePosition[yj] + leftOffset;
				 xi < nX - rightOffset; xi++) {
				// We are only interested in the helium near the surface
				if ((grid[xi] + grid[xi + 1]) / 2.0 -
						grid[surfacePosition[yj] + 1] >
					2.0)
					continue;

				// Check if we are on the right processor
				if (xi >= localXS && xi < localXS + localXM &&
					yj >= bottomOffset && yj < nY - topOffset) {
					// Get the concentrations at this grid point
					auto concOffset =
						subview(concs, yj, xi, Kokkos::ALL).view();

					// Sum the total atom concentration
					atomConc += psiNetwork.getTotalTrappedHeliumConcentration(
									concOffset, 0) *
						(grid[xi + 1] - grid[xi]);
				}
			}

			// Share the concentration with all the processes
			totalAtomConc = 0.0;
			auto xolotlComm = util::getMPIComm();
			MPI_Allreduce(
				&atomConc, &totalAtomConc, 1, MPI_DOUBLE, MPI_SUM, xolotlComm);

			// Set the disappearing rate in the modified TM handler
			psiNetwork.updateTrapMutationDisappearingRate(totalAtomConc);
		}

		// Set the grid position
		gridPosition[1] = yj * hY;

		// Initialize the advection and temperature handlers which depend
		// on the surface position at Y
		advectionHandlers[0]->setLocation(
			grid[surfacePosition[yj] + 1] - grid[1]);

		for (auto xi = localXS; xi < localXS + localXM; xi++) {
			// Boundary conditions
			// Everything to the left of the surface is empty
			if (xi < surfacePosition[yj] + leftOffset ||
				xi > nX - 1 - rightOffset || yj < bottomOffset ||
				yj > nY - 1 - topOffset) {
				valIndex += 5 * nDiff;
				valIndex += 2 * nAdvec * advectionHandlers.size();
				valIndex += nNetworkEntries;
				continue;
			}
			// Free surface GB
			if (std::find_if(begin(gbVector), end(gbVector), [=](auto&& pair) {
					return xi == pair[0] && yj == pair[1];
				}) != end(gbVector)) {
				valIndex += 5 * nDiff;
				valIndex += 2 * nAdvec * advectionHandlers.size();
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
			diffusionHandler->computePartialsForDiffusion(network,
				subview(vals, std::make_pair(valIndex, valIndex + 5 * nDiff)),
				hxLeft, hxRight, xi - localXS, sy, yj - localYS);
			valIndex += 5 * nDiff;

			// Get the partial derivatives for the advection
			// Set the grid position
			gridPosition[0] = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];
			for (auto l = 0; l < advectionHandlers.size(); l++) {
				advectionHandlers[l]->computePartialsForAdvection(network,
					subview(
						vals, std::make_pair(valIndex, valIndex + 2 * nAdvec)),
					gridPosition, hxLeft, hxRight, xi - localXS, hY,
					yj - localYS);
				valIndex += 2 * nAdvec;
			}

			// Get the concentations
			auto concOffset = subview(concs, yj, xi, Kokkos::ALL).view();

			// ----- Take care of the reactions for all the reactants -----

			auto surfacePos = grid[surfacePosition[yj] + 1];
			auto curXPos = (grid[xi] + grid[xi + 1]) / 2.0;
			auto prevXPos = (grid[xi - 1] + grid[xi]) / 2.0;
			auto curDepth = curXPos - surfacePos;
			auto curSpacing = curXPos - prevXPos;

			// Compute all the partial derivatives for the reactions
			partialDerivativeCounter->increment();
			partialDerivativeTimer->start();
			network.computeAllPartials(concOffset,
				subview(
					vals, std::make_pair(valIndex, valIndex + nNetworkEntries)),
				xi + 1 - localXS, curDepth, curSpacing);
			partialDerivativeTimer->stop();
			valIndex += nNetworkEntries;
		}
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
} /* end namespace solver */
} /* end namespace xolotl */
