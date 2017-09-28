// Includes
#include <PetscSolver3DHandler.h>
#include <MathUtils.h>
#include <Constants.h>

namespace xolotlSolver {

void PetscSolver3DHandler::createSolverContext(DM &da) {
	PetscErrorCode ierr;

	// Initialize the all reactants pointer
	allReactants = network->getAll();

	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network->getDOF();

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	ierr = DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_MIRROR,
			DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_STAR, nX,
			nY, nZ, PETSC_DECIDE,
			PETSC_DECIDE, PETSC_DECIDE, dof, 1, NULL, NULL, NULL, &da);
	checkPetscError(ierr, "PetscSolver3DHandler::createSolverContext: "
			"DMDACreate3d failed.");
	ierr = DMSetFromOptions(da);
	checkPetscError(ierr,
			"PetscSolver3DHandler::createSolverContext: DMSetFromOptions failed.");
	ierr = DMSetUp(da);
	checkPetscError(ierr,
			"PetscSolver3DHandler::createSolverContext: DMSetUp failed.");

	// Set the position of the surface
	// Loop on Y
	for (int j = 0; j < nY; j++) {
		// Create a one dimensional vector to store the surface indices
		// for a given Y position
		std::vector<int> tempPosition;

		// Loop on Z
		for (int k = 0; k < nZ; k++) {
			tempPosition.push_back(0);
			if (movingSurface)
				tempPosition[k] = (int) (nX * portion / 100.0);
		}

		// Add tempPosition to the surfacePosition
		surfacePosition.push_back(tempPosition);
	}

	// Generate the grid in the x direction
	generateGrid(nX, hX, surfacePosition[0][0]);

	// Now that the grid was generated, we can update the surface position
	// if we are using a restart file
	int tempTimeStep = -2;
	bool hasConcentrations = false;
	if (!networkName.empty())
		hasConcentrations = xolotlCore::HDF5Utils::hasConcentrationGroup(
			networkName, tempTimeStep);

	// Get the actual surface position if concentrations were stored
	if (hasConcentrations) {
		auto surfaceIndices = xolotlCore::HDF5Utils::readSurface3D(networkName,
				tempTimeStep);

		// Set the actual surface positions
		for (int i = 0; i < surfaceIndices.size(); i++) {
			for (int j = 0; j < surfaceIndices[0].size(); j++) {
				surfacePosition[i][j] = surfaceIndices[i][j];
			}
		}
	}

	// Initialize the surface of the first advection handler corresponding to the
	// advection toward the surface (or a dummy one if it is deactivated)
	advectionHandlers[0]->setLocation(grid[surfacePosition[0][0]]);

	// Set the size of the partial derivatives vector
	reactingPartialsForCluster.resize(dof, 0.0);

	/*  The only spatial coupling in the Jacobian is due to diffusion.
	 *  The ofill (thought of as a dof by dof 2d (row-oriented) array represents
	 *  the nonzero coupling between degrees of freedom at one point with degrees
	 *  of freedom on the adjacent point to the left or right. A 1 at i,j in the
	 *  ofill array indicates that the degree of freedom i at a point is coupled
	 *  to degree of freedom j at the adjacent point.
	 *  In this case ofill has only a few diagonal entries since the only spatial
	 *  coupling is regular diffusion.
	 */
	PetscInt *ofill, *dfill;
	ierr = PetscMalloc(dof * dof * sizeof(PetscInt), &ofill);
	checkPetscError(ierr, "PetscSolver3DHandler::createSolverContext: "
			"PetscMalloc (ofill) failed.");
	ierr = PetscMalloc(dof * dof * sizeof(PetscInt), &dfill);
	checkPetscError(ierr, "PetscSolver3DHandler::createSolverContext: "
			"PetscMalloc (dfill) failed.");
	ierr = PetscMemzero(ofill, dof * dof * sizeof(PetscInt));
	checkPetscError(ierr, "PetscSolver3DHandler::createSolverContext: "
			"PetscMemzero (ofill) failed.");
	ierr = PetscMemzero(dfill, dof * dof * sizeof(PetscInt));
	checkPetscError(ierr, "PetscSolver3DHandler::createSolverContext: "
			"PetscMemzero (dfill) failed.");

	// Initialize the temperature handler
	temperatureHandler->initializeTemperature(network, ofill, dfill);

	// Fill ofill, the matrix of "off-diagonal" elements that represents diffusion
	diffusionHandler->initializeOFill(network, ofill);
	// Loop on the advection handlers to account the other "off-diagonal" elements
	for (int i = 0; i < advectionHandlers.size(); i++) {
		advectionHandlers[i]->initialize(network, ofill);
	}

	// Initialize the modified trap-mutation handler because it adds connectivity
	mutationHandler->initialize(network, grid, nY, hY, nZ, hZ);
	mutationHandler->initializeIndex3D(surfacePosition, network,
			advectionHandlers, grid, nY, hY, nZ, hZ);

	// Get the diagonal fill
	network->getDiagonalFill(dfill);

	// Load up the block fills
	ierr = DMDASetBlockFills(da, dfill, ofill);
	checkPetscError(ierr, "PetscSolver3DHandler::createSolverContext: "
			"DMDASetBlockFills failed.");

	// Free the temporary fill arrays
	ierr = PetscFree(ofill);
	checkPetscError(ierr, "PetscSolver3DHandler::createSolverContext: "
			"PetscFree (ofill) failed.");
	ierr = PetscFree(dfill);
	checkPetscError(ierr, "PetscSolver3DHandler::createSolverContext: "
			"PetscFree (dfill) failed.");

	// Initialize the arrays for the reaction partial derivatives
	reactionVals = new PetscScalar[dof * dof];
	reactionIndices = new PetscInt[dof * dof];

	return;
}

void PetscSolver3DHandler::initializeConcentration(DM &da, Vec &C) {
	PetscErrorCode ierr;

	// Pointer for the concentration vector
	PetscScalar ****concentrations = nullptr;
	ierr = DMDAVecGetArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolver3DHandler::initializeConcentration: "
			"DMDAVecGetArrayDOF failed.");

	// Get the local boundaries
	PetscInt xs, xm, ys, ym, zs, zm;
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
	checkPetscError(ierr, "PetscSolver3DHandler::initializeConcentration: "
			"DMDAGetCorners failed.");

	// Get the last time step written in the HDF5 file
	int tempTimeStep = -2;
	bool hasConcentrations = false;
	if (!networkName.empty())
		hasConcentrations = xolotlCore::HDF5Utils::hasConcentrationGroup(
			networkName, tempTimeStep);

	// Give the surface position to the temperature handler
	temperatureHandler->updateSurfacePosition(grid[surfacePosition[0][0]]);

	// Initialize the flux handler
	fluxHandler->initializeFluxHandler(network, surfacePosition[0][0], grid);

	// Initialize the grid for the diffusion
	diffusionHandler->initializeDiffusionGrid(advectionHandlers, grid, nY, hY,
			nZ, hZ);

	// Initialize the grid for the advection
	advectionHandlers[0]->initializeAdvectionGrid(advectionHandlers, grid, nY,
			hY, nZ, hZ);

	// Pointer for the concentration vector at a specific grid point
	PetscScalar *concOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network->getDOF();

	// Get the single vacancy ID
	auto singleVacancyCluster = network->get(xolotlCore::vType, 1);
	int vacancyIndex = -1;
	if (singleVacancyCluster)
		vacancyIndex = singleVacancyCluster->getId() - 1;

	// Loop on all the grid points
	for (PetscInt k = zs; k < zs + zm; k++) {
		for (PetscInt j = ys; j < ys + ym; j++) {
			for (PetscInt i = xs; i < xs + xm; i++) {
				concOffset = concentrations[k][j][i];

				// Loop on all the clusters to initialize at 0.0
				for (int n = 0; n < dof - 1; n++) {
					concOffset[n] = 0.0;
				}

				// Temperature
				std::vector<double> gridPosition = { grid[i], 0.0, 0.0 };
				concOffset[dof - 1] = temperatureHandler->getTemperature(gridPosition,
						0.0);

				// Initialize the vacancy concentration
				if (i > surfacePosition[j][k]
						&& vacancyIndex > 0 && !hasConcentrations) {
					concOffset[vacancyIndex] = initialVConc;
				}
			}
		}
	}

	// If the concentration must be set from the HDF5 file
	if (hasConcentrations) {
		// Loop on the full grid
		for (PetscInt k = 0; k < nZ; k++) {
			for (PetscInt j = 0; j < nY; j++) {
				for (PetscInt i = 0; i < nX; i++) {
					// Read the concentrations from the HDF5 file
					auto concVector = xolotlCore::HDF5Utils::readGridPoint(
							networkName, tempTimeStep, i, j, k);

					// Change the concentration only if we are on the locally
					// owned part of the grid
					if (i >= xs && i < xs + xm && j >= ys && j < ys + ym
							&& k >= zs && k < zs + zm) {
						concOffset = concentrations[k][j][i];
						// Loop on the concVector size
						for (unsigned int l = 0; l < concVector.size(); l++) {
							concOffset[(int) concVector.at(l).at(0)] =
									concVector.at(l).at(1);
						}
					}
				}
			}
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolver3DHandler::initializeConcentration: "
			"DMDAVecRestoreArrayDOF failed.");

	return;
}

void PetscSolver3DHandler::updateConcentration(TS &ts, Vec &localC, Vec &F,
		PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the local data vector from PETSc
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver3DHandler::updateConcentration: "
			"TSGetDM failed.");

	// Pointers to the PETSc arrays that start at the beginning (xs, ys, zs) of the
	// local array!
	PetscScalar ****concs = nullptr, ****updatedConcs = nullptr;
	// Get pointers to vector data
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver3DHandler::updateConcentration: "
			"DMDAVecGetArrayDOFRead (localC) failed.");
	ierr = DMDAVecGetArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolver3DHandler::updateConcentration: "
			"DMDAVecGetArrayDOF (F) failed.");

	// Get local grid boundaries
	PetscInt xs, xm, ys, ym, zs, zm;
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
	checkPetscError(ierr, "PetscSolver3DHandler::updateConcentration: "
			"DMDAGetCorners failed.");

	// The following pointers are set to the first position in the conc or
	// updatedConc arrays that correspond to the beginning of the data for the
	// current grid point. They are accessed just like regular arrays.
	PetscScalar *concOffset = nullptr, *updatedConcOffset = nullptr;

	// Set some step size variable
	double sy = 1.0 / (hY * hY);
	double sz = 1.0 / (hZ * hZ);

	// Declarations for variables used in the loop
	double **concVector = new double*[7];
	std::vector<double> gridPosition = { 0.0, 0.0, 0.0 }, incidentFluxVector;
	double atomConc = 0.0, totalAtomConc = 0.0;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network->getDOF();

	// Loop over grid points
	for (PetscInt zk = 0; zk < nZ; zk++) {
		for (PetscInt yj = 0; yj < nY; yj++) {

			// Compute the total concentration of atoms contained in bubbles
			atomConc = 0.0;

			// Loop over grid points
			for (int xi = surfacePosition[yj][zk] + 1; xi < nX; xi++) {
				// We are only interested in the helium near the surface
				if (grid[xi] - grid[surfacePosition[yj][zk]] > 2.0)
					continue;

				// Check if we are on the right processor
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym
						&& zk >= zs && zk < zs + zm) {
					// Get the concentrations at this grid point
					concOffset = concs[zk][yj][xi];
					// Copy data into the PSIClusterReactionNetwork
					network->updateConcentrationsFromArray(concOffset);

					// Sum the total atom concentration
					atomConc += network->getTotalTrappedAtomConcentration()
							* (grid[xi] - grid[xi - 1]);
				}
			}

			// Share the concentration with all the processes
			totalAtomConc = 0.0;
			MPI_Allreduce(&atomConc, &totalAtomConc, 1, MPI_DOUBLE, MPI_SUM,
					MPI_COMM_WORLD);

			// Set the disappearing rate in the modified TM handler
			mutationHandler->updateDisappearingRate(totalAtomConc);

			// Skip if we are not on the right process
			if (yj < ys || yj >= ys + ym || zk < zs || zk >= zs + zm)
				continue;

			// Set the grid position
			gridPosition[1] = yj * hY;
			gridPosition[2] = zk * hZ;

			// Initialize the flux, advection, and temperature handlers which depend
			// on the surface position at Y
			fluxHandler->initializeFluxHandler(network, surfacePosition[yj][zk],
					grid);
			advectionHandlers[0]->setLocation(grid[surfacePosition[yj][zk]]);
			temperatureHandler->updateSurfacePosition(grid[surfacePosition[yj][zk]]);

			for (PetscInt xi = xs; xi < xs + xm; xi++) {
				// Compute the old and new array offsets
				concOffset = concs[zk][yj][xi];
				updatedConcOffset = updatedConcs[zk][yj][xi];

				// Boundary conditions
				// Everything to the left of the surface is empty
				if (xi <= surfacePosition[yj][zk]) {
					for (int i = 0; i < dof; i++) {
						updatedConcOffset[i] = 1.0 * concOffset[i];
					}

					continue;
				}

				// Set the grid position
				gridPosition[0] = grid[xi];

				// Fill the concVector with the pointer to the middle, left,
				// right, bottom, top, front, and back grid points
				concVector[0] = concOffset; // middle
				concVector[1] = concs[zk][yj][xi - 1]; // left
				concVector[2] = concs[zk][yj][xi + 1]; // right
				concVector[3] = concs[zk][yj - 1][xi]; // bottom
				concVector[4] = concs[zk][yj + 1][xi]; // top
				concVector[5] = concs[zk - 1][yj][xi]; // front
				concVector[6] = concs[zk + 1][yj][xi]; // back

				// Get the temperature from the temperature handler
				temperatureHandler->setTemperature(concOffset);
				double temperature = temperatureHandler->getTemperature(gridPosition,
						ftime);

				// Update the network if the temperature changed
				if (!xolotlCore::equal(temperature, lastTemperature)) {
					network->setTemperature(temperature);
					// Update the modified trap-mutation rate that depends on the
					// network reaction rates
					mutationHandler->updateTrapMutationRate(network);
					lastTemperature = temperature;
				}

				// Copy data into the ReactionNetwork so that it can
				// compute the fluxes properly. The network is only used to compute the
				// fluxes and hold the state data from the last time step. I'm reusing
				// it because it cuts down on memory significantly (about 400MB per
				// grid point) at the expense of being a little tricky to comprehend.
				network->updateConcentrationsFromArray(concOffset);

				// ----- Account for flux of incoming particles -----
				fluxHandler->computeIncidentFlux(ftime, updatedConcOffset, xi, surfacePosition[yj][zk]);

				// ---- Compute the temperature over the locally owned part of the grid -----
				temperatureHandler->computeTemperature(concVector, updatedConcOffset,
						grid[xi] - grid[xi - 1], grid[xi + 1] - grid[xi]);

				// ---- Compute diffusion over the locally owned part of the grid -----
				diffusionHandler->computeDiffusion(network, concVector,
						updatedConcOffset, grid[xi] - grid[xi - 1],
						grid[xi + 1] - grid[xi], xi, sy, yj, sz, zk);

				// ---- Compute advection over the locally owned part of the grid -----
				for (int i = 0; i < advectionHandlers.size(); i++) {
					advectionHandlers[i]->computeAdvection(network,
							gridPosition, concVector, updatedConcOffset,
							grid[xi] - grid[xi - 1], grid[xi + 1] - grid[xi],
							xi, hY, yj, hZ, zk);
				}

				// ----- Compute the modified trap-mutation over the locally owned part of the grid -----
				mutationHandler->computeTrapMutation(network, concOffset,
						updatedConcOffset, xi, yj, zk);

				// ----- Compute the reaction fluxes over the locally owned part of the grid -----
				network->computeAllFluxes(updatedConcOffset);
			}
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver3DHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOFRead (localC) failed.");
	ierr = DMDAVecRestoreArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolver3DHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOF (F) failed.");
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr, "PetscSolver3DHandler::updateConcentration: "
			"DMRestoreLocalVector failed.");

	// Clear memory
	delete[] concVector;

	return;
}

void PetscSolver3DHandler::computeOffDiagonalJacobian(TS &ts, Vec &localC,
		Mat &J, PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver3DHandler::computeOffDiagonalJacobian: "
			"TSGetDM failed.");

	// Setup some step size variables
	double sy = 1.0 / (hY * hY);
	double sz = 1.0 / (hZ * hZ);

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network->getDOF();

	// Pointers to the PETSc arrays that start at the beginning (xs) of the
	// local array!
	PetscScalar ****concs = nullptr;
	// Get pointers to vector data
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver3DHandler::computeOffDiagonalJacobian: "
			"DMDAVecGetArrayDOFRead (localC) failed.");

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset = nullptr;

	// Get local grid boundaries
	PetscInt xs, xm, ys, ym, zs, zm;
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
	checkPetscError(ierr, "PetscSolver3DHandler::computeOffDiagonalJacobian: "
			"DMDAGetCorners failed.");

	// Get the total number of diffusing clusters
	const int nDiff = max(diffusionHandler->getNumberOfDiffusing(), 1);

	// Get the total number of advecting clusters
	int nAdvec = 0;
	for (int l = 0; l < advectionHandlers.size(); l++) {
		int n = advectionHandlers[l]->getNumberOfAdvecting();
		if (n > nAdvec)
			nAdvec = n;
	}

	// Arguments for MatSetValuesStencil called below
	MatStencil row, cols[7];
	PetscScalar diffVals[7 * nDiff];
	PetscInt diffIndices[nDiff];
	PetscScalar advecVals[2 * nAdvec];
	PetscInt advecIndices[nAdvec];
	std::vector<double> gridPosition = { 0.0, 0.0, 0.0 };

	/*
	 Loop over grid points computing Jacobian terms for diffusion and advection
	 at each grid point
	 */
	for (PetscInt zk = zs; zk < zs + zm; zk++) {
		// Set the grid position
		gridPosition[2] = zk * hZ;
		for (PetscInt yj = ys; yj < ys + ym; yj++) {
			// Set the grid position
			gridPosition[1] = yj * hY;

			// Initialize the advection and temperature handlers which depend
			// on the surface position at Y
			advectionHandlers[0]->setLocation(grid[surfacePosition[yj][zk]]);
			temperatureHandler->updateSurfacePosition(grid[surfacePosition[yj][zk]]);

			for (PetscInt xi = xs; xi < xs + xm; xi++) {
				// Boundary conditions
				// Everything to the left of the surface is empty
				if (xi <= surfacePosition[yj][zk])
					continue;

				// Set the grid position
				gridPosition[0] = grid[xi];

				// Get the temperature from the temperature handler
				concOffset = concs[zk][yj][xi];
				temperatureHandler->setTemperature(concOffset);
				double temperature = temperatureHandler->getTemperature(gridPosition,
						ftime);

				// Update the network if the temperature changed
				if (!xolotlCore::equal(temperature, lastTemperature)) {
					network->setTemperature(temperature);
					// Update the modified trap-mutation rate that depends on the
					// network reaction rates
					mutationHandler->updateTrapMutationRate(network);
					lastTemperature = temperature;
				}

				// Get the partial derivatives for the temperature
				temperatureHandler->computePartialsForTemperature(diffVals, diffIndices,
						grid[xi] - grid[xi - 1], grid[xi + 1] - grid[xi]);

				// Set grid coordinate and component number for the row
				row.i = xi;
				row.j = yj;
				row.k = zk;
				row.c = diffIndices[0];

				// Set grid coordinates and component numbers for the columns
				// corresponding to the middle, left, and right grid points
				cols[0].i = xi; // middle
				cols[0].j = yj;
				cols[0].k = zk;
				cols[0].c = diffIndices[0];
				cols[1].i = xi - 1; // left
				cols[1].j = yj;
				cols[1].k = zk;
				cols[1].c = diffIndices[0];
				cols[2].i = xi + 1; // right
				cols[2].j = yj;
				cols[2].k = zk;
				cols[2].c = diffIndices[0];

				ierr = MatSetValuesStencil(J, 1, &row, 3, cols, diffVals, ADD_VALUES);
				checkPetscError(ierr,
						"PetscSolver3DHandler::computeOffDiagonalJacobian: "
								"MatSetValuesStencil (temperature) failed.");

				// Get the partial derivatives for the diffusion
				diffusionHandler->computePartialsForDiffusion(network, diffVals,
						diffIndices, grid[xi] - grid[xi - 1],
						grid[xi + 1] - grid[xi], xi, sy, yj, sz, zk);

				// Loop on the number of diffusion cluster to set the values in the Jacobian
				for (int i = 0; i < nDiff; i++) {
					// Set grid coordinate and component number for the row
					row.i = xi;
					row.j = yj;
					row.k = zk;
					row.c = diffIndices[i];

					// Set grid coordinates and component numbers for the columns
					// corresponding to the middle, left, right, bottom, top, front,
					// and back grid points
					cols[0].i = xi; // middle
					cols[0].j = yj;
					cols[0].k = zk;
					cols[0].c = diffIndices[i];
					cols[1].i = xi - 1; // left
					cols[1].j = yj;
					cols[1].k = zk;
					cols[1].c = diffIndices[i];
					cols[2].i = xi + 1; // right
					cols[2].j = yj;
					cols[2].k = zk;
					cols[2].c = diffIndices[i];
					cols[3].i = xi; // bottom
					cols[3].j = yj - 1;
					cols[3].k = zk;
					cols[3].c = diffIndices[i];
					cols[4].i = xi; // top
					cols[4].j = yj + 1;
					cols[4].k = zk;
					cols[4].c = diffIndices[i];
					cols[5].i = xi; // front
					cols[5].j = yj;
					cols[5].k = zk - 1;
					cols[5].c = diffIndices[i];
					cols[6].i = xi; // back
					cols[6].j = yj;
					cols[6].k = zk + 1;
					cols[6].c = diffIndices[i];

					ierr = MatSetValuesStencil(J, 1, &row, 7, cols,
							diffVals + (7 * i), ADD_VALUES);
					checkPetscError(ierr,
							"PetscSolver3DHandler::computeOffDiagonalJacobian: "
									"MatSetValuesStencil (diffusion) failed.");
				}

				// Get the partial derivatives for the advection
				for (int l = 0; l < advectionHandlers.size(); l++) {
					advectionHandlers[l]->computePartialsForAdvection(network,
							advecVals, advecIndices, gridPosition,
							grid[xi] - grid[xi - 1], grid[xi + 1] - grid[xi],
							xi, hY, yj, hZ, zk);

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
						row.k = zk;
						row.c = advecIndices[i];

						// If we are on the sink, the partial derivatives are not the same
						// Both sides are giving their concentrations to the center
						if (advectionHandlers[l]->isPointOnSink(gridPosition)) {
							cols[0].i = xi - advecStencil[0]; // left?
							cols[0].j = yj - advecStencil[1]; // bottom?
							cols[0].k = zk - advecStencil[2]; // back?
							cols[0].c = advecIndices[i];
							cols[1].i = xi + advecStencil[0]; // right?
							cols[1].j = yj + advecStencil[1]; // top?
							cols[1].k = zk + advecStencil[2]; // front?
							cols[1].c = advecIndices[i];
						} else {
							// Set grid coordinates and component numbers for the columns
							// corresponding to the middle and other grid points
							cols[0].i = xi; // middle
							cols[0].j = yj;
							cols[0].k = zk;
							cols[0].c = advecIndices[i];
							cols[1].i = xi + advecStencil[0]; // left or right?
							cols[1].j = yj + advecStencil[1]; // bottom or top?
							cols[1].k = zk + advecStencil[2]; // back or front?
							cols[1].c = advecIndices[i];
						}

						// Update the matrix
						ierr = MatSetValuesStencil(J, 1, &row, 2, cols,
								advecVals + (2 * i), ADD_VALUES);
						checkPetscError(ierr,
								"PetscSolver3DHandler::computeOffDiagonalJacobian: "
										"MatSetValuesStencil (advection) failed.");
					}
				}
			}
		}
	}

	return;
}

void PetscSolver3DHandler::computeDiagonalJacobian(TS &ts, Vec &localC, Mat &J,
		PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver3DHandler::computeDiagonalJacobian: "
			"TSGetDM failed.");

	// Get pointers to vector data
	PetscScalar ****concs = nullptr;
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver3DHandler::computeDiagonalJacobian: "
			"DMDAVecGetArrayDOFRead failed.");

	// Get local grid boundaries
	PetscInt xs, xm, ys, ym, zs, zm;
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
	checkPetscError(ierr, "PetscSolver3DHandler::computeDiagonalJacobian: "
			"DMDAGetCorners failed.");

	// The degree of freedom is the size of the network
	const int dof = network->getDOF();

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset = nullptr;

	// Arguments for MatSetValuesStencil called below
	MatStencil rowId;
	MatStencil colIds[dof];
	int pdColIdsVectorSize = 0;
	PetscInt reactionSize[dof];

	// Store the total number of He clusters in the network for the
	// modified trap-mutation
	int nHelium = network->getAll(xolotlCore::heType).size();

	// Declarations for variables used in the loop
	double atomConc = 0.0, totalAtomConc = 0.0;
	std::vector<double> gridPosition = { 0.0, 0.0, 0.0 };

	// Loop over the grid points
	for (PetscInt zk = 0; zk < nZ; zk++) {
		for (PetscInt yj = 0; yj < nY; yj++) {

			// Compute the total concentration of atoms contained in bubbles
			atomConc = 0.0;

			// Loop over grid points
			for (int xi = surfacePosition[yj][zk] + 1; xi < nX; xi++) {
				// We are only interested in the helium near the surface
				if (grid[xi] - grid[surfacePosition[yj][zk]] > 2.0)
					continue;

				// Check if we are on the right processor
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym
						&& zk >= zs && zk < zs + zm) {
					// Get the concentrations at this grid point
					concOffset = concs[zk][yj][xi];
					// Copy data into the PSIClusterReactionNetwork
					network->updateConcentrationsFromArray(concOffset);

					// Sum the total atom concentration
					atomConc += network->getTotalTrappedAtomConcentration()
							* (grid[xi] - grid[xi - 1]);
				}
			}

			// Share the concentration with all the processes
			totalAtomConc = 0.0;
			MPI_Allreduce(&atomConc, &totalAtomConc, 1, MPI_DOUBLE, MPI_SUM,
					MPI_COMM_WORLD);

			// Set the disappearing rate in the modified TM handler
			mutationHandler->updateDisappearingRate(totalAtomConc);

			// Skip if we are not on the right process
			if (yj < ys || yj >= ys + ym || zk < zs || zk >= zs + zm)
				continue;

			// Set the grid position
			gridPosition[1] = yj * hY;
			gridPosition[2] = zk * hZ;

			for (PetscInt xi = xs; xi < xs + xm; xi++) {
				// Boundary conditions
				// Everything to the left of the surface is empty
				if (xi <= surfacePosition[yj][zk])
					continue;

				// Set the grid position
				gridPosition[0] = grid[xi];

				// Get the temperature from the temperature handler
				concOffset = concs[zk][yj][xi];
				temperatureHandler->setTemperature(concOffset);
				double temperature = temperatureHandler->getTemperature(gridPosition,
						ftime);

				// Update the network if the temperature changed
				if (!xolotlCore::equal(temperature, lastTemperature)) {
					network->setTemperature(temperature);
					// Update the modified trap-mutation rate that depends on the
					// network reaction rates
					mutationHandler->updateTrapMutationRate(network);
					lastTemperature = temperature;
				}

				// Copy data into the ReactionNetwork so that it can
				// compute the new concentrations.
				network->updateConcentrationsFromArray(concOffset);

				// ----- Take care of the reactions for all the reactants -----

				// Compute all the partial derivatives for the reactions
				network->computeAllPartials(reactionVals, reactionIndices,
						reactionSize);

				// Update the column in the Jacobian that represents each DOF
				for (int i = 0; i < dof - 1; i++) {
					// Set grid coordinate and component number for the row
					rowId.i = xi;
					rowId.j = yj;
					rowId.k = zk;
					rowId.c = i;

					// Number of partial derivatives
					pdColIdsVectorSize = reactionSize[i];
					// Loop over the list of column ids
					for (int j = 0; j < pdColIdsVectorSize; j++) {
						// Set grid coordinate and component number for a column in the list
						colIds[j].i = xi;
						colIds[j].j = yj;
						colIds[j].k = zk;
						colIds[j].c = reactionIndices[i * dof + j];
						// Get the partial derivative from the array of all of the partials
						reactingPartialsForCluster[j] =
								reactionVals[i * dof + j];
					}
					// Update the matrix
					ierr = MatSetValuesStencil(J, 1, &rowId, pdColIdsVectorSize,
							colIds, reactingPartialsForCluster.data(),
							ADD_VALUES);
					checkPetscError(ierr,
							"PetscSolver1DHandler::computeDiagonalJacobian: "
									"MatSetValuesStencil (reactions) failed.");
				}

				// ----- Take care of the modified trap-mutation for all the reactants -----

				// Arguments for MatSetValuesStencil called below
				MatStencil row, col;
				PetscScalar mutationVals[3 * nHelium];
				PetscInt mutationIndices[3 * nHelium];

				// Compute the partial derivative from modified trap-mutation at this grid point
				int nMutating = mutationHandler->computePartialsForTrapMutation(
						network, mutationVals, mutationIndices, xi, yj, zk);

				// Loop on the number of helium undergoing trap-mutation to set the values
				// in the Jacobian
				for (int i = 0; i < nMutating; i++) {
					// Set grid coordinate and component number for the row and column
					// corresponding to the helium cluster
					row.i = xi;
					row.j = yj;
					row.k = zk;
					row.c = mutationIndices[3 * i];
					col.i = xi;
					col.j = yj;
					col.k = zk;
					col.c = mutationIndices[3 * i];

					ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
							mutationVals + (3 * i), ADD_VALUES);
					checkPetscError(ierr,
							"PetscSolver3DHandler::computeDiagonalJacobian: "
									"MatSetValuesStencil (He trap-mutation) failed.");

					// Set component number for the row
					// corresponding to the HeV cluster created through trap-mutation
					row.c = mutationIndices[(3 * i) + 1];

					ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
							mutationVals + (3 * i) + 1, ADD_VALUES);
					checkPetscError(ierr,
							"PetscSolver3DHandler::computeDiagonalJacobian: "
									"MatSetValuesStencil (HeV trap-mutation) failed.");

					// Set component number for the row
					// corresponding to the interstitial created through trap-mutation
					row.c = mutationIndices[(3 * i) + 2];

					ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
							mutationVals + (3 * i) + 2, ADD_VALUES);
					checkPetscError(ierr,
							"PetscSolver3DHandler::computeDiagonalJacobian: "
									"MatSetValuesStencil (I trap-mutation) failed.");
				}
			}
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver3DHandler::computeDiagonalJacobian: "
			"DMDAVecRestoreArrayDOFRead failed.");
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr, "PetscSolver3DHandler::computeDiagonalJacobian: "
			"DMRestoreLocalVector failed.");

	return;
}

} /* end namespace xolotlSolver */
