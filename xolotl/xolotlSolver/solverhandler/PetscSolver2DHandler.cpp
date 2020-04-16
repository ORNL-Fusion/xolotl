// Includes
#include <PetscSolver2DHandler.h>
#include <MathUtils.h>
#include <Constants.h>
#include <experimental/PSIReactionNetwork.h>

namespace xolotlSolver {

void PetscSolver2DHandler::createSolverContext(DM &da) {
	PetscErrorCode ierr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = expNetwork.getDOF();

	// Set the position of the surface
	for (int j = 0; j < nY; j++) {
		surfacePosition.push_back(0);
		if (movingSurface)
			surfacePosition[j] = (int) (nX * portion / 100.0);
	}

	// Generate the grid in the x direction
	generateGrid(nX, hX, surfacePosition[0]);

	// Now that the grid was generated, we can update the surface position
	// if we are using a restart file
	if (not networkName.empty() and movingSurface) {

		xolotlCore::XFile xfile(networkName);
		auto concGroup =
				xfile.getGroup<xolotlCore::XFile::ConcentrationGroup>();
		if (concGroup and concGroup->hasTimesteps()) {

			auto tsGroup = concGroup->getLastTimestepGroup();
			assert(tsGroup);
			auto surfaceIndices = tsGroup->readSurface2D();

			// Set the actual surface positions
			for (int i = 0; i < surfaceIndices.size(); i++) {
				surfacePosition[i] = surfaceIndices[i];
			}
		}
	}

	// Prints the grid on one process
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);
	if (procId == 0) {
		for (int i = 1; i < grid.size() - 1; i++) {
			std::cout << grid[i] - grid[surfacePosition[0] + 1] << " ";
		}
		std::cout << std::endl;
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_MIRROR,
			DM_BOUNDARY_PERIODIC, DMDA_STENCIL_STAR, nX, nY, PETSC_DECIDE,
			PETSC_DECIDE, dof + 1, 1, NULL, NULL, &da);
	checkPetscError(ierr, "PetscSolver2DHandler::createSolverContext: "
			"DMDACreate2d failed.");
	ierr = DMSetFromOptions(da);
	checkPetscError(ierr,
			"PetscSolver2DHandler::createSolverContext: DMSetFromOptions failed.");
	ierr = DMSetUp(da);
	checkPetscError(ierr,
			"PetscSolver2DHandler::createSolverContext: DMSetUp failed.");

	// Initialize the surface of the first advection handler corresponding to the
	// advection toward the surface (or a dummy one if it is deactivated)
	advectionHandlers[0]->setLocation(grid[surfacePosition[0] + 1] - grid[1]);

	/*  The only spatial coupling in the Jacobian is due to diffusion.
	 *  The ofill (thought of as a dof by dof 2d (row-oriented) array represents
	 *  the nonzero coupling between degrees of freedom at one point with degrees
	 *  of freedom on the adjacent point to the left or right. A 1 at i,j in the
	 *  ofill array indicates that the degree of freedom i at a point is coupled
	 *  to degree of freedom j at the adjacent point.
	 *  In this case ofill has only a few diagonal entries since the only spatial
	 *  coupling is regular diffusion.
	 */
	xolotlCore::IReactionNetwork::SparseFillMap ofill;

	// Initialize the temperature handler
	temperatureHandler->initializeTemperature(dof, ofill, dfill);

	// Fill ofill, the matrix of "off-diagonal" elements that represents diffusion
	diffusionHandler->initializeOFill(expNetwork, ofill);
	// Loop on the advection handlers to account the other "off-diagonal" elements
	for (int i = 0; i < advectionHandlers.size(); i++) {
		advectionHandlers[i]->initialize(expNetwork, ofill);
	}

	// Get the local boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr, "PetscSolver2DHandler::initializeConcentration: "
			"DMDAGetCorners failed.");

	// Initialize the modified trap-mutation handler because it adds connectivity
	mutationHandler->initialize(expNetwork, dfill, xm, ym);
	mutationHandler->initializeIndex2D(surfacePosition, expNetwork,
			advectionHandlers, grid, xm, xs, ym, hY, ys);

	// Tell the network the number of grid points on this process with ghosts
	// TODO: do we need the ghost points?
	expNetwork.setGridSize(xm + 2);

	// Get the diagonal fill
	auto nPartials = expNetwork.getDiagonalFill(dfill);

	// Load up the block fills
	auto dfillsparse = ConvertToPetscSparseFillMap(dof + 1, dfill);
	auto ofillsparse = ConvertToPetscSparseFillMap(dof + 1, ofill);
	ierr = DMDASetBlockFillsSparse(da, dfillsparse.data(), ofillsparse.data());
	checkPetscError(ierr, "PetscSolver2DHandler::createSolverContext: "
			"DMDASetBlockFills failed.");

	// Initialize the arrays for the reaction partial derivatives
	expVals = Kokkos::View<double*>("solverPartials", nPartials);

	// Set the size of the partial derivatives vectors
	reactingPartialsForCluster.resize(dof, 0.0);

	return;
}

void PetscSolver2DHandler::initializeConcentration(DM &da, Vec &C) {
	PetscErrorCode ierr;

	// Pointer for the concentration vector
	PetscScalar ***concentrations = nullptr;
	ierr = DMDAVecGetArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolver2DHandler::initializeConcentration: "
			"DMDAVecGetArrayDOF failed.");

	// Get the local boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr, "PetscSolver2DHandler::initializeConcentration: "
			"DMDAGetCorners failed.");

	// Initialize the last temperature at each grid point on this process
	for (int i = 0; i < xm + 2; i++) {
		temperature.push_back(0.0);
	}

	// Get the last time step written in the HDF5 file
	bool hasConcentrations = false;
	std::unique_ptr<xolotlCore::XFile> xfile;
	std::unique_ptr<xolotlCore::XFile::ConcentrationGroup> concGroup;
	if (not networkName.empty()) {
		xfile.reset(new xolotlCore::XFile(networkName));
		concGroup = xfile->getGroup<xolotlCore::XFile::ConcentrationGroup>();
		hasConcentrations = (concGroup and concGroup->hasTimesteps());
	}

	// Give the surface position to the temperature handler
	temperatureHandler->updateSurfacePosition(surfacePosition[0]);

	// Initialize the flux handler
	fluxHandler->initializeFluxHandler(expNetwork, surfacePosition[0], grid);

	// Initialize the grid for the diffusion
	diffusionHandler->initializeDiffusionGrid(advectionHandlers, grid, xm, xs,
			ym, hY, ys);

	// Initialize the grid for the advection
	advectionHandlers[0]->initializeAdvectionGrid(advectionHandlers, grid, xm,
			xs, ym, hY, ys);

	// Pointer for the concentration vector at a specific grid point
	PetscScalar *concOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	// + the super clusters
	const int dof = expNetwork.getDOF();

	// Get the single vacancy ID
	auto singleVacancyCluster = expNetwork.getSingleVacancy();
	auto vacancyIndex = NetworkType::invalidIndex();
	if (singleVacancyCluster.getId() != NetworkType::invalidIndex())
		vacancyIndex = singleVacancyCluster.getId();

	// Loop on all the grid points
	for (PetscInt j = ys; j < ys + ym; j++)
		for (PetscInt i = xs - 1; i <= xs + xm; i++) {
			// Temperature
			xolotlCore::Point<3> gridPosition { 0.0, j * hY, 0.0 };
			if (i < 0)
				gridPosition[0] =
						(grid[0] - grid[surfacePosition[j] + 1])
								/ (grid[grid.size() - 1]
										- grid[surfacePosition[j] + 1]);
			else
				gridPosition[0] =
						((grid[i] + grid[i + 1]) / 2.0
								- grid[surfacePosition[j] + 1])
								/ (grid[grid.size() - 1]
										- grid[surfacePosition[j] + 1]);
			auto temp = temperatureHandler->getTemperature(gridPosition, 0.0);
			temperature[i - xs + 1] = temp;

			// Boundary conditions
			if (i < xs || i >= xs + xm)
				continue;

			concOffset = concentrations[j][i];
			concOffset[dof] = temp;

			// Loop on all the clusters to initialize at 0.0
			for (int n = 0; n < dof; n++) {
				concOffset[n] = 0.0;
			}

			// Initialize the vacancy concentration
			if (i >= surfacePosition[j] + leftOffset
					and vacancyIndex != NetworkType::invalidIndex()
					and not hasConcentrations and i < nX - rightOffset
					and j >= bottomOffset and j < nY - topOffset) {
				concOffset[vacancyIndex] = initialVConc;
			}
		}

	// If the concentration must be set from the HDF5 file
	if (hasConcentrations) {

		assert(concGroup);
		auto tsGroup = concGroup->getLastTimestepGroup();
		assert(tsGroup);

		// Loop on the full grid
		for (PetscInt j = 0; j < nY; j++)
			for (PetscInt i = 0; i < nX; i++) {
				// Read the concentrations from the HDF5 file
				auto concVector = tsGroup->readGridPoint(i, j);

				// Change the concentration only if we are on the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym) {
					concOffset = concentrations[j][i];
					// Loop on the concVector size
					for (unsigned int l = 0; l < concVector.size(); l++) {
						concOffset[(int) concVector.at(l).at(0)] =
								concVector.at(l).at(1);
					}
					// Get the temperature
					double temp = concVector.at(concVector.size() - 1).at(1);
					temperature[i - xs + 1] = temp;
				}
			}
	}

	// Update the network with the temperature
	expNetwork.setTemperatures(temperature);
	expNetwork.syncClusterDataOnHost();
	// Update the modified trap-mutation rate
	// that depends on the network reaction rates
	mutationHandler->updateTrapMutationRate(expNetwork.getLargestRate());

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolver2DHandler::initializeConcentration: "
			"DMDAVecRestoreArrayDOF failed.");

	return;
}

void PetscSolver2DHandler::updateConcentration(TS &ts, Vec &localC, Vec &F,
		PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the local data vector from PETSc
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr,
			"PetscSolver2DHandler::updateConcentration: TSGetDM failed.");

	// Pointers to the PETSc arrays that start at the beginning (xs, ys) of the
	// local array!
	PetscScalar ***concs = nullptr, ***updatedConcs = nullptr;
	// Get pointers to vector data
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: "
			"DMDAVecGetArrayDOFRead (localC) failed.");
	ierr = DMDAVecGetArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: "
			"DMDAVecGetArrayDOF (F) failed.");

	// Get local grid boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: "
			"DMDAGetCorners failed.");

	// The following pointers are set to the first position in the conc or
	// updatedConc arrays that correspond to the beginning of the data for the
	// current grid point. They are accessed just like regular arrays.
	PetscScalar *concOffset = nullptr, *updatedConcOffset = nullptr;

	// Set some step size variable
	double sy = 1.0 / (hY * hY);

	// Declarations for variables used in the loop
	double **concVector = new double*[5];
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };
	std::vector<double> incidentFluxVector;
	double atomConc = 0.0, totalAtomConc = 0.0;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = expNetwork.getDOF();

	// Loop over grid points first for the temperature, including the ghost points in X
	for (PetscInt yj = ys; yj < ys + ym; yj++) {
		temperatureHandler->updateSurfacePosition(surfacePosition[yj]);
		bool tempHasChanged = false;
		for (PetscInt xi = xs - 1; xi <= xs + xm; xi++) {

			// Heat condition
			if (xi == surfacePosition[yj] && xi >= xs && xi < xs + xm) {
				// Compute the old and new array offsets
				concOffset = concs[yj][xi];
				updatedConcOffset = updatedConcs[yj][xi];

				// Fill the concVector with the pointer to the middle, left, and right grid points
				concVector[0] = concOffset; // middle
				concVector[1] = concs[yj][xi - 1]; // left
				concVector[2] = concs[yj][xi + 1]; // right
				concVector[3] = concs[yj - 1][xi]; // bottom
				concVector[4] = concs[yj + 1][xi]; // top

				// Compute the left and right hx
				double hxLeft = 0.0, hxRight = 0.0;
				if (xi - 1 >= 0 && xi < nX) {
					hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
					hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
				} else if (xi == -1) {
					hxLeft = 0.0;
					hxRight = (grid[xi + 2] + grid[xi + 1]) / 2.0;
				} else if (xi - 1 < 0) {
					hxLeft = (grid[xi + 1] + grid[xi]) / 2.0;
					hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
				} else {
					hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
					hxRight = (grid[xi + 1] - grid[xi]) / 2;
				}

				temperatureHandler->computeTemperature(concVector,
						updatedConcOffset, hxLeft, hxRight, xi, sy, yj);
			}

			// Boundary conditions
			// Everything to the left of the surface is empty
			if (xi < surfacePosition[yj] + leftOffset
					|| xi > nX - 1 - rightOffset) {
				continue;
			}
			// Free surface GB
			bool skip = false;
			for (auto &pair : gbVector) {
				if (xi == std::get<0>(pair) && yj == std::get<1>(pair)) {
					skip = true;
					break;
				}
			}
			if (skip)
				continue;
			// Compute the old and new array offsets
			concOffset = concs[yj][xi];
			updatedConcOffset = updatedConcs[yj][xi];

			// Fill the concVector with the pointer to the middle, left, and right grid points
			concVector[0] = concOffset; // middle
			concVector[1] = concs[yj][xi - 1]; // left
			concVector[2] = concs[yj][xi + 1]; // right
			concVector[3] = concs[yj - 1][xi]; // bottom
			concVector[4] = concs[yj + 1][xi]; // top

			// Compute the left and right hx
			double hxLeft = 0.0, hxRight = 0.0;
			if (xi - 1 >= 0 && xi < nX) {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			} else if (xi - 1 < 0) {
				hxLeft = (grid[xi + 1] + grid[xi]) / 2.0;
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			} else {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = (grid[xi + 1] - grid[xi]) / 2;
			}

			// Set the grid fraction
			gridPosition[0] = ((grid[xi] + grid[xi + 1]) / 2.0
					- grid[surfacePosition[yj] + 1])
					/ (grid[grid.size() - 1] - grid[surfacePosition[yj] + 1]);
			gridPosition[1] = yj / nY;

			// Get the temperature from the temperature handler
			temperatureHandler->setTemperature(concOffset);
			double temp = temperatureHandler->getTemperature(gridPosition,
					ftime);

			// Update the network if the temperature changed
			if (std::fabs(temperature[xi + 1 - xs] - temp) > 0.1) {
				temperature[xi + 1 - xs] = temp;
				tempHasChanged = true;
			}

			// ---- Compute the temperature over the locally owned part of the grid -----
			if (xi >= xs && xi < xs + xm) {
				temperatureHandler->computeTemperature(concVector,
						updatedConcOffset, hxLeft, hxRight, xi, sy, yj);
			}
		}

		// TODO: it is updated T more than once per MPI process in preparation
		// of T depending on more than X
		if (tempHasChanged) {
			// Update the network with the temperature
			expNetwork.setTemperatures(temperature);
			expNetwork.syncClusterDataOnHost();
			// Update the modified trap-mutation rate
			// that depends on the network reaction rates
			// TODO: is this just the local largest rate? Is it correct?
			mutationHandler->updateTrapMutationRate(
					expNetwork.getLargestRate());
		}
	}

	// Loop over grid points
	for (PetscInt yj = bottomOffset; yj < nY - topOffset; yj++) {

		// Computing the trapped atom concentration is only needed for the attenuation
		if (useAttenuation) {
			// Compute the total concentration of atoms contained in bubbles
			atomConc = 0.0;

			// Loop over grid points
			for (int xi = surfacePosition[yj] + leftOffset;
					xi < nX - rightOffset; xi++) {
				// We are only interested in the helium near the surface
				if ((grid[xi] + grid[xi + 1]) / 2.0
						- grid[surfacePosition[yj] + 1] > 2.0)
					continue;

				// Check if we are on the right processor
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
					// Get the concentrations at this grid point
					concOffset = concs[yj][xi];

					// Sum the total atom concentration
					using NetworkType =
					xolotlCore::experimental::PSIReactionNetwork<xolotlCore::experimental::PSIFullSpeciesList>;
					using Spec = typename NetworkType::Species;
					using HostUnmanaged =
					Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
					auto hConcs = HostUnmanaged(concOffset, dof);
					auto dConcs = Kokkos::View<double*>("Concentrations", dof);
					deep_copy(dConcs, hConcs);
					// TODO: how to not have to cast the network here?
					auto &psiNetwork = dynamic_cast<NetworkType&>(expNetwork);
					atomConc += psiNetwork.getTotalTrappedAtomConcentration(
							dConcs, Spec::He, 0) * (grid[xi + 1] - grid[xi]);
				}
			}

			// Share the concentration with all the processes
			totalAtomConc = 0.0;
			MPI_Allreduce(&atomConc, &totalAtomConc, 1, MPI_DOUBLE, MPI_SUM,
					MPI_COMM_WORLD);

			// Set the disappearing rate in the modified TM handler
			mutationHandler->updateDisappearingRate(totalAtomConc);
		}

		// Skip if we are not on the right process
		if (yj < ys || yj >= ys + ym)
			continue;

		// Set the grid position
		gridPosition[1] = yj * hY;

		// Initialize the flux, advection handlers which depend
		// on the surface position at Y
		fluxHandler->initializeFluxHandler(expNetwork, surfacePosition[yj],
				grid);
		advectionHandlers[0]->setLocation(
				grid[surfacePosition[yj] + 1] - grid[1]);

		for (PetscInt xi = xs; xi < xs + xm; xi++) {
			// Compute the old and new array offsets
			concOffset = concs[yj][xi];
			updatedConcOffset = updatedConcs[yj][xi];

			// Fill the concVector with the pointer to the middle, left, right, bottom, and top grid points
			concVector[0] = concOffset; // middle
			concVector[1] = concs[yj][xi - 1]; // left
			concVector[2] = concs[yj][xi + 1]; // right
			concVector[3] = concs[yj - 1][xi]; // bottom
			concVector[4] = concs[yj + 1][xi]; // top

			// Compute the left and right hx
			double hxLeft = 0.0, hxRight = 0.0;
			if (xi - 1 >= 0 && xi < nX) {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			} else if (xi - 1 < 0) {
				hxLeft = (grid[xi + 1] + grid[xi]) / 2.0;
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			} else {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = (grid[xi + 1] - grid[xi]) / 2;
			}

			// Boundary conditions
			// Everything to the left of the surface is empty
			if (xi < surfacePosition[yj] + leftOffset
					|| xi > nX - 1 - rightOffset || yj < bottomOffset
					|| yj > nY - 1 - topOffset) {
				continue;
			}
			// Free surface GB
			bool skip = false;
			for (auto &pair : gbVector) {
				if (xi == std::get<0>(pair) && yj == std::get<1>(pair)) {
					skip = true;
					break;
				}
			}
			if (skip)
				continue;

			// ----- Account for flux of incoming particles -----
			fluxHandler->computeIncidentFlux(ftime, updatedConcOffset, xi,
					surfacePosition[yj]);

			// ---- Compute diffusion over the locally owned part of the grid -----
			diffusionHandler->computeDiffusion(expNetwork, concVector,
					updatedConcOffset, hxLeft, hxRight, xi - xs, sy, yj - ys);

			// ---- Compute advection over the locally owned part of the grid -----
			// Set the grid position
			gridPosition[0] = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];
			for (int i = 0; i < advectionHandlers.size(); i++) {
				advectionHandlers[i]->computeAdvection(expNetwork, gridPosition,
						concVector, updatedConcOffset, hxLeft, hxRight, xi - xs,
						hY, yj - ys);
			}

			// ----- Compute the modified trap-mutation over the locally owned part of the grid -----
			mutationHandler->computeTrapMutation(expNetwork, concOffset,
					updatedConcOffset, xi - xs, yj - ys);

			// ----- Compute the reaction fluxes over the locally owned part of the grid -----
			using HostUnmanaged =
			Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
			auto hConcs = HostUnmanaged(concOffset, dof + 1);
			auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
			deep_copy(dConcs, hConcs);
			auto hFlux = HostUnmanaged(updatedConcOffset, dof + 1);
			auto dFlux = Kokkos::View<double*>("Fluxes", dof + 1);
			deep_copy(dFlux, hFlux);
			fluxCounter->increment();
			fluxTimer->start();
			expNetwork.computeAllFluxes(dConcs, dFlux, xi + 1 - xs);
			fluxTimer->stop();
			deep_copy(hFlux, dFlux);
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOFRead (localC) failed.");
	ierr = DMDAVecRestoreArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOF (F) failed.");

	// Clear memory
	delete[] concVector;

	return;
}

void PetscSolver2DHandler::computeOffDiagonalJacobian(TS &ts, Vec &localC,
		Mat &J, PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver2DHandler::computeOffDiagonalJacobian: "
			"TSGetDM failed.");

	// Setup some step size variables
	double sy = 1.0 / (hY * hY);

	// Degrees of freedom is the total number of clusters in the network
	const int dof = expNetwork.getDOF();

	// Pointers to the PETSc arrays that start at the beginning (xs) of the
	// local array!
	PetscScalar ***concs = nullptr;
	// Get pointers to vector data
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver2DHandler::computeOffDiagonalJacobian: "
			"DMDAVecGetArrayDOFRead (localC) failed.");

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset = nullptr;

	// Get local grid boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr, "PetscSolver2DHandler::computeOffDiagonalJacobian: "
			"DMDAGetCorners failed.");

	// Get the total number of diffusing clusters
	const int nDiff = max(diffusionHandler->getNumberOfDiffusing(), 0);

	// Get the total number of advecting clusters
	int nAdvec = 0;
	for (int l = 0; l < advectionHandlers.size(); l++) {
		int n = advectionHandlers[l]->getNumberOfAdvecting();
		if (n > nAdvec)
			nAdvec = n;
	}

	// Arguments for MatSetValuesStencil called below
	MatStencil row, cols[5];
	PetscScalar tempVals[5];
	PetscInt tempIndices[1];
	PetscScalar diffVals[5 * nDiff];
	PetscInt diffIndices[nDiff];
	PetscScalar advecVals[2 * nAdvec];
	PetscInt advecIndices[nAdvec];
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };

	/*
	 Loop over grid points for the temperature, including ghosts
	 */
	for (PetscInt yj = ys; yj < ys + ym; yj++) {
		temperatureHandler->updateSurfacePosition(surfacePosition[yj]);
		bool tempHasChanged = false;
		for (PetscInt xi = xs - 1; xi <= xs + xm; xi++) {
			// Compute the left and right hx
			double hxLeft = 0.0, hxRight = 0.0;
			if (xi - 1 >= 0 && xi < nX) {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			} else if (xi == -1) {
				hxLeft = 0.0;
				hxRight = (grid[xi + 2] + grid[xi + 1]) / 2.0;
			} else if (xi - 1 < 0) {
				hxLeft = (grid[xi + 1] + grid[xi]) / 2.0;
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			} else {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = (grid[xi + 1] - grid[xi]) / 2;
			}

			// Heat condition
			if (xi == surfacePosition[yj] && xi >= xs && xi < xs + xm) {
				// Get the partial derivatives for the temperature
				auto setValues =
						temperatureHandler->computePartialsForTemperature(
								tempVals, tempIndices, hxLeft, hxRight, xi, sy,
								yj);

				if (setValues) {
					// Set grid coordinate and component number for the row
					row.i = xi;
					row.j = yj;
					row.c = tempIndices[0];

					// Set grid coordinates and component numbers for the columns
					// corresponding to the middle, left, and right grid points
					cols[0].i = xi; // middle
					cols[0].j = yj;
					cols[0].c = tempIndices[0];
					cols[1].i = xi - 1; // left
					cols[1].j = yj;
					cols[1].c = tempIndices[0];
					cols[2].i = xi + 1; // right
					cols[2].j = yj;
					cols[2].c = tempIndices[0];
					cols[3].i = xi; // bottom
					cols[3].j = yj - 1;
					cols[3].c = tempIndices[0];
					cols[4].i = xi; // top
					cols[4].j = yj + 1;
					cols[4].c = tempIndices[0];

					ierr = MatSetValuesStencil(J, 1, &row, 5, cols, tempVals,
							ADD_VALUES);
					checkPetscError(ierr,
							"PetscSolver2DHandler::computeOffDiagonalJacobian: "
									"MatSetValuesStencil (temperature) failed.");
				}
			}

			// Boundary conditions
			// Everything to the left of the surface is empty
			if (xi < surfacePosition[yj] + leftOffset
					|| xi > nX - 1 - rightOffset)
				continue;
			// Free surface GB
			bool skip = false;
			for (auto &pair : gbVector) {
				if (xi == std::get<0>(pair) && yj == std::get<1>(pair)) {
					skip = true;
					break;
				}
			}
			if (skip)
				continue;

			// Get the concentrations at this grid point
			concOffset = concs[yj][xi];

			// Set the grid fraction
			gridPosition[0] = ((grid[xi] + grid[xi + 1]) / 2.0
					- grid[surfacePosition[yj] + 1])
					/ (grid[grid.size() - 1] - grid[surfacePosition[yj] + 1]);
			gridPosition[1] = yj / nY;

			// Get the temperature from the temperature handler
			temperatureHandler->setTemperature(concOffset);
			double temp = temperatureHandler->getTemperature(gridPosition,
					ftime);

			// Update the network if the temperature changed
			if (std::fabs(temperature[xi + 1 - xs] - temp) > 0.1) {
				temperature[xi + 1 - xs] = temp;
				tempHasChanged = true;
			}

			// Get the partial derivatives for the temperature
			if (xi >= xs && xi < xs + xm) {
				auto setValues =
						temperatureHandler->computePartialsForTemperature(
								tempVals, tempIndices, hxLeft, hxRight, xi, sy,
								yj);

				if (setValues) {
					// Set grid coordinate and component number for the row
					row.i = xi;
					row.j = yj;
					row.c = tempIndices[0];

					// Set grid coordinates and component numbers for the columns
					// corresponding to the middle, left, and right grid points
					cols[0].i = xi; // middle
					cols[0].j = yj;
					cols[0].c = tempIndices[0];
					cols[1].i = xi - 1; // left
					cols[1].j = yj;
					cols[1].c = tempIndices[0];
					cols[2].i = xi + 1; // right
					cols[2].j = yj;
					cols[2].c = tempIndices[0];
					cols[3].i = xi; // bottom
					cols[3].j = yj - 1;
					cols[3].c = tempIndices[0];
					cols[4].i = xi; // top
					cols[4].j = yj + 1;
					cols[4].c = tempIndices[0];

					ierr = MatSetValuesStencil(J, 1, &row, 5, cols, tempVals,
							ADD_VALUES);
					checkPetscError(ierr,
							"PetscSolver2DHandler::computeOffDiagonalJacobian: "
									"MatSetValuesStencil (temperature) failed.");
				}
			}
		}

		if (tempHasChanged) {
			// Update the network with the temperature
			expNetwork.setTemperatures(temperature);
			expNetwork.syncClusterDataOnHost();
		}
	}

	/*
	 Loop over grid points computing Jacobian terms for diffusion and advection
	 at each grid point
	 */
	for (PetscInt yj = ys; yj < ys + ym; yj++) {
		// Set the grid position
		gridPosition[1] = yj * hY;

		// Initialize the advection and temperature handlers which depend
		// on the surface position at Y
		advectionHandlers[0]->setLocation(
				grid[surfacePosition[yj] + 1] - grid[1]);

		for (PetscInt xi = xs; xi < xs + xm; xi++) {
			// Compute the left and right hx
			double hxLeft = 0.0, hxRight = 0.0;
			if (xi - 1 >= 0 && xi < nX) {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			} else if (xi - 1 < 0) {
				hxLeft = (grid[xi + 1] + grid[xi]) / 2.0;
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			} else {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = (grid[xi + 1] - grid[xi]) / 2;
			}

			// Boundary conditions
			// Everything to the left of the surface is empty
			if (xi < surfacePosition[yj] + leftOffset
					|| xi > nX - 1 - rightOffset || yj < bottomOffset
					|| yj > nY - 1 - topOffset)
				continue;
			// Free surface GB
			bool skip = false;
			for (auto &pair : gbVector) {
				if (xi == std::get<0>(pair) && yj == std::get<1>(pair)) {
					skip = true;
					break;
				}
			}
			if (skip)
				continue;

			// Get the partial derivatives for the diffusion
			diffusionHandler->computePartialsForDiffusion(expNetwork, diffVals,
					diffIndices, hxLeft, hxRight, xi - xs, sy, yj - ys);

			// Loop on the number of diffusion cluster to set the values in the Jacobian
			for (int i = 0; i < nDiff; i++) {
				// Set grid coordinate and component number for the row
				row.i = xi;
				row.j = yj;
				row.c = diffIndices[i];

				// Set grid coordinates and component numbers for the columns
				// corresponding to the middle, left, right, bottom, and top grid points
				cols[0].i = xi; // middle
				cols[0].j = yj;
				cols[0].c = diffIndices[i];
				cols[1].i = xi - 1; // left
				cols[1].j = yj;
				cols[1].c = diffIndices[i];
				cols[2].i = xi + 1; // right
				cols[2].j = yj;
				cols[2].c = diffIndices[i];
				cols[3].i = xi; // bottom
				cols[3].j = yj - 1;
				cols[3].c = diffIndices[i];
				cols[4].i = xi; // top
				cols[4].j = yj + 1;
				cols[4].c = diffIndices[i];

				ierr = MatSetValuesStencil(J, 1, &row, 5, cols,
						diffVals + (5 * i), ADD_VALUES);
				checkPetscError(ierr,
						"PetscSolver2DHandler::computeOffDiagonalJacobian: "
								"MatSetValuesStencil (diffusion) failed.");
			}

			// Get the partial derivatives for the advection
			// Set the grid position
			gridPosition[0] = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];
			for (int l = 0; l < advectionHandlers.size(); l++) {
				advectionHandlers[l]->computePartialsForAdvection(expNetwork,
						advecVals, advecIndices, gridPosition, hxLeft, hxRight,
						xi - xs, hY, yj - ys);

				// Get the stencil indices to know where to put the partial derivatives in the Jacobian
				auto advecStencil =
						advectionHandlers[l]->getStencilForAdvection(
								gridPosition);

				// Get the number of advecting clusters
				nAdvec = advectionHandlers[l]->getNumberOfAdvecting();

				// Loop on the number of advecting cluster to set the values in the Jacobian
				for (int i = 0; i < nAdvec; i++) {
					// Set grid coordinate and component number for the row
					row.i = xi;
					row.j = yj;
					row.c = advecIndices[i];

					// If we are on the sink, the partial derivatives are not the same
					// Both sides are giving their concentrations to the center
					if (advectionHandlers[l]->isPointOnSink(gridPosition)) {
						cols[0].i = xi - advecStencil[0]; // left?
						cols[0].j = yj - advecStencil[1]; // bottom?
						cols[0].c = advecIndices[i];
						cols[1].i = xi + advecStencil[0]; // right?
						cols[1].j = yj + advecStencil[1]; // top?
						cols[1].c = advecIndices[i];
					} else {
						// Set grid coordinates and component numbers for the columns
						// corresponding to the middle and other grid points
						cols[0].i = xi; // middle
						cols[0].j = yj;
						cols[0].c = advecIndices[i];
						cols[1].i = xi + advecStencil[0]; // left or right?
						cols[1].j = yj + advecStencil[1]; // bottom or top?
						cols[1].c = advecIndices[i];
					}

					// Update the matrix
					ierr = MatSetValuesStencil(J, 1, &row, 2, cols,
							advecVals + (2 * i), ADD_VALUES);
					checkPetscError(ierr,
							"PetscSolver2DHandler::computeOffDiagonalJacobian: "
									"MatSetValuesStencil (advection) failed.");
				}
			}
		}
	}

	// Restore the array
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver2DHandler::computeOffDiagonalJacobian: "
			"DMDAVecRestoreArrayDOFRead (localC) failed.");

	return;
}

void PetscSolver2DHandler::computeDiagonalJacobian(TS &ts, Vec &localC, Mat &J,
		PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
			"TSGetDM failed.");

	// Get pointers to vector data
	PetscScalar ***concs = nullptr;
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
			"DMDAVecGetArrayDOFRead failed.");

	// Get local grid boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
			"DMDAGetCorners failed.");

	// The degree of freedom is the size of the network
	const int dof = expNetwork.getDOF();

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset = nullptr;

	// Arguments for MatSetValuesStencil called below
	MatStencil rowId;
	MatStencil colIds[dof];
	int pdColIdsVectorSize = 0;

	// Declarations for variables used in the loop
	double atomConc = 0.0, totalAtomConc = 0.0;
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };

	// Loop over grid points first for the temperature, including the ghost points in X
	for (PetscInt yj = ys; yj < ys + ym; yj++) {
		temperatureHandler->updateSurfacePosition(surfacePosition[yj]);
		bool tempHasChanged = false;
		for (PetscInt xi = xs - 1; xi <= xs + xm; xi++) {
			// Boundary conditions
			// Everything to the left of the surface is empty
			if (xi < surfacePosition[yj] + leftOffset
					|| xi > nX - 1 - rightOffset) {
				continue;
			}
			// Free surface GB
			bool skip = false;
			for (auto &pair : gbVector) {
				if (xi == std::get<0>(pair) && yj == std::get<1>(pair)) {
					skip = true;
					break;
				}
			}
			if (skip)
				continue;

			// Get the concentrations at this grid point
			concOffset = concs[yj][xi];

			// Set the grid fraction
			gridPosition[0] = ((grid[xi] + grid[xi + 1]) / 2.0
					- grid[surfacePosition[yj] + 1])
					/ (grid[grid.size() - 1] - grid[surfacePosition[yj] + 1]);
			gridPosition[1] = yj / nY;

			// Get the temperature from the temperature handler
			temperatureHandler->setTemperature(concOffset);
			double temp = temperatureHandler->getTemperature(gridPosition,
					ftime);

			// Update the network if the temperature changed
			if (std::fabs(temperature[xi + 1 - xs] - temp) > 0.1) {
				temperature[xi + 1 - xs] = temp;
				tempHasChanged = true;
			}
		}

		// TODO: it is updated T more than once per MPI process in preparation
		// of T depending on more than X
		if (tempHasChanged) {
			// Update the network with the temperature
			expNetwork.setTemperatures(temperature);
			expNetwork.syncClusterDataOnHost();
			// Update the modified trap-mutation rate
			// that depends on the network reaction rates
			// TODO: is this just the local largest rate? Is it correct?
			mutationHandler->updateTrapMutationRate(
					expNetwork.getLargestRate());
		}
	}

	// Loop over the grid points
	for (PetscInt yj = bottomOffset; yj < nY - topOffset; yj++) {

		// Computing the trapped atom concentration is only needed for the attenuation
		if (useAttenuation) {
			// Compute the total concentration of atoms contained in bubbles
			atomConc = 0.0;

			// Loop over grid points
			for (int xi = surfacePosition[yj] + leftOffset;
					xi < nX - rightOffset; xi++) {
				// We are only interested in the helium near the surface
				if ((grid[xi] + grid[xi + 1]) / 2.0
						- grid[surfacePosition[yj] + 1] > 2.0)
					continue;

				// Check if we are on the right processor
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
					// Get the concentrations at this grid point
					concOffset = concs[yj][xi];

					// Sum the total atom concentration
					using NetworkType =
					xolotlCore::experimental::PSIReactionNetwork<xolotlCore::experimental::PSIFullSpeciesList>;
					using Spec = typename NetworkType::Species;
					using HostUnmanaged =
					Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
					auto hConcs = HostUnmanaged(concOffset, dof);
					auto dConcs = Kokkos::View<double*>("Concentrations", dof);
					deep_copy(dConcs, hConcs);
					// TODO: how to not have to cast the network here?
					auto &psiNetwork = dynamic_cast<NetworkType&>(expNetwork);
					atomConc += psiNetwork.getTotalTrappedAtomConcentration(
							dConcs, Spec::He, 0) * (grid[xi + 1] - grid[xi]);
				}
			}

			// Share the concentration with all the processes
			totalAtomConc = 0.0;
			MPI_Allreduce(&atomConc, &totalAtomConc, 1, MPI_DOUBLE, MPI_SUM,
					MPI_COMM_WORLD);

			// Set the disappearing rate in the modified TM handler
			mutationHandler->updateDisappearingRate(totalAtomConc);
		}

		// Skip if we are not on the right process
		if (yj < ys || yj >= ys + ym)
			continue;

		// Set the grid position
		gridPosition[1] = yj * hY;

		for (PetscInt xi = xs; xi < xs + xm; xi++) {
			// Boundary conditions
			// Everything to the left of the surface is empty
			if (xi < surfacePosition[yj] + leftOffset
					|| xi > nX - 1 - rightOffset || yj < bottomOffset
					|| yj > nY - 1 - topOffset)
				continue;
			// Free surface GB
			bool skip = false;
			for (auto &pair : gbVector) {
				if (xi == std::get<0>(pair) && yj == std::get<1>(pair)) {
					skip = true;
					break;
				}
			}
			if (skip)
				continue;

			// Set the grid fraction
			gridPosition[0] = ((grid[xi] + grid[xi + 1]) / 2.0
					- grid[surfacePosition[yj] + 1])
					/ (grid[grid.size() - 1] - grid[surfacePosition[yj] + 1]);

			// Get the concentations
			concOffset = concs[yj][xi];

			// ----- Take care of the reactions for all the reactants -----

			// Compute all the partial derivatives for the reactions
			using HostUnmanaged =
			Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
			auto hConcs = HostUnmanaged(concOffset, dof + 1);
			auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
			deep_copy(dConcs, hConcs);
			partialDerivativeCounter->increment();
			partialDerivativeTimer->start();
			expNetwork.computeAllPartials(dConcs, expVals, xi + 1 - xs);
			partialDerivativeTimer->stop();
			auto hPartials = create_mirror_view(expVals);
			deep_copy(hPartials, expVals);

			// Variable for the loop on reactants
			int startingIdx = 0;
			// Update the column in the Jacobian that represents each DOF
			for (int i = 0; i < dof; i++) {
				// Set grid coordinate and component number for the row
				rowId.i = xi;
				rowId.j = yj;
				rowId.c = i;

				// Number of partial derivatives
				auto rowIter = dfill.find(i);
				if (rowIter != dfill.end()) {
					const auto &row = rowIter->second;
					pdColIdsVectorSize = row.size();

					// Loop over the list of column ids
					for (int j = 0; j < pdColIdsVectorSize; j++) {
						// Set grid coordinate and component number for a column in the list
						colIds[j].i = xi;
						colIds[j].j = yj;
						colIds[j].c = row[j];
						// Get the partial derivative from the array of all of the partials
						reactingPartialsForCluster[j] = hPartials(
								startingIdx + j);
					}
					// Update the matrix
					ierr = MatSetValuesStencil(J, 1, &rowId, pdColIdsVectorSize,
							colIds, reactingPartialsForCluster.data(),
							ADD_VALUES);
					checkPetscError(ierr,
							"PetscSolver2DHandler::computeDiagonalJacobian: "
									"MatSetValuesStencil (reactions) failed.");

					// Increase the starting index
					startingIdx += pdColIdsVectorSize;
				}
			}

			// ----- Take care of the modified trap-mutation for all the reactants -----

			// Store the total number of He clusters in the network for the
			// modified trap-mutation
			int nHelium = mutationHandler->getNumberOfMutating();

			// Arguments for MatSetValuesStencil called below
			MatStencil row, col;
			PetscScalar mutationVals[3 * nHelium];
			PetscInt mutationIndices[3 * nHelium];

			// Compute the partial derivative from modified trap-mutation at this grid point
			int nMutating = mutationHandler->computePartialsForTrapMutation(
					expNetwork, concOffset, mutationVals, mutationIndices,
					xi - xs, yj - ys);

			// Loop on the number of helium undergoing trap-mutation to set the values
			// in the Jacobian
			for (int i = 0; i < nMutating; i++) {
				// Set grid coordinate and component number for the row and column
				// corresponding to the helium cluster
				row.i = xi;
				row.j = yj;
				row.c = mutationIndices[3 * i];
				col.i = xi;
				col.j = yj;
				col.c = mutationIndices[3 * i];

				ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
						mutationVals + (3 * i), ADD_VALUES);
				checkPetscError(ierr,
						"PetscSolver2DHandler::computeDiagonalJacobian: "
								"MatSetValuesStencil (He trap-mutation) failed.");

				// Set component number for the row
				// corresponding to the HeV cluster created through trap-mutation
				row.c = mutationIndices[(3 * i) + 1];

				ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
						mutationVals + (3 * i) + 1, ADD_VALUES);
				checkPetscError(ierr,
						"PetscSolver2DHandler::computeDiagonalJacobian: "
								"MatSetValuesStencil (HeV trap-mutation) failed.");

				// Set component number for the row
				// corresponding to the interstitial created through trap-mutation
				row.c = mutationIndices[(3 * i) + 2];

				ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
						mutationVals + (3 * i) + 2, ADD_VALUES);
				checkPetscError(ierr,
						"PetscSolver2DHandler::computeDiagonalJacobian: "
								"MatSetValuesStencil (I trap-mutation) failed.");
			}
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
			"DMDAVecRestoreArrayDOFRead failed.");

	return;
}

} /* end namespace xolotlSolver */
