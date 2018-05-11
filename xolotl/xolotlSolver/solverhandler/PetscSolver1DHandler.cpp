// Includes
#include <PetscSolver1DHandler.h>
#include <MathUtils.h>
#include <Constants.h>

namespace xolotlSolver {

void PetscSolver1DHandler::createSolverContext(DM &da) {

	PetscErrorCode ierr;
	// Recompute Ids and network size and redefine the connectivities
	network.reinitializeConnectivities();

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network.getDOF();

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_MIRROR, nX, dof, 1,
	NULL, &da);
	checkPetscError(ierr, "PetscSolver1DHandler::createSolverContext: "
			"DMDACreate1d failed.");
	ierr = DMSetFromOptions(da);
	checkPetscError(ierr,
			"PetscSolver1DHandler::createSolverContext: DMSetFromOptions failed.");
	ierr = DMSetUp(da);
	checkPetscError(ierr,
			"PetscSolver1DHandler::createSolverContext: DMSetUp failed.");

	// Set the position of the surface
	surfacePosition = 0;
	if (movingSurface)
		surfacePosition = (int) (nX * portion / 100.0);

	// Generate the grid in the x direction
	generateGrid(nX, hX, surfacePosition);

	// Now that the grid was generated, we can update the surface position
	// if we are using a restart file
	int tempTimeStep = -2;
	bool hasConcentrations = false;
	if (!networkName.empty())
		hasConcentrations = xolotlCore::HDF5Utils::hasConcentrationGroup(
				networkName, tempTimeStep);
	if (hasConcentrations) {
		surfacePosition = xolotlCore::HDF5Utils::readSurface1D(networkName,
				tempTimeStep);
	}

	// Initialize the surface of the first advection handler corresponding to the
	// advection toward the surface (or a dummy one if it is deactivated)
	advectionHandlers[0]->setLocation(grid[surfacePosition + 1] - grid[1]);

	// Prints the grid on one process
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);
	if (procId == 0) {
		for (int i = 0; i < grid.size() - 1; i++) {
			std::cout << grid[i + 1] - grid[surfacePosition + 1] << " ";
		}
		std::cout << std::endl;
	}

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
	xolotlCore::IReactionNetwork::SparseFillMap ofill;
	xolotlCore::IReactionNetwork::SparseFillMap dfill;

	// Initialize the temperature handler
	temperatureHandler->initializeTemperature(network, ofill, dfill);

	// Fill ofill, the matrix of "off-diagonal" elements that represents diffusion
	diffusionHandler->initializeOFill(network, ofill);
	// Loop on the advection handlers to account the other "off-diagonal" elements
	for (int i = 0; i < advectionHandlers.size(); i++) {
		advectionHandlers[i]->initialize(network, ofill);
	}

	// Initialize the modified trap-mutation handler here
	// because it adds connectivity
	mutationHandler->initialize(network, grid);
	mutationHandler->initializeIndex1D(surfacePosition, network,
			advectionHandlers, grid);

	// Get the diagonal fill
	network.getDiagonalFill(dfill);

	// Load up the block fills
	auto ofillsparse = ConvertToPetscSparseFillMap(dof, ofill);
	auto dfillsparse = ConvertToPetscSparseFillMap(dof, dfill);
	ierr = DMDASetBlockFillsSparse(da, dfillsparse.data(), ofillsparse.data());
	checkPetscError(ierr, "PetscSolver1DHandler::createSolverContext: "
			"DMDASetBlockFills failed.");

	// Initialize the arrays for the reaction partial derivatives
	reactionSize.resize(dof);
	reactionStartingIdx.resize(dof);
	auto nPartials = network.initPartialsSizes(reactionSize,
			reactionStartingIdx);

	reactionIndices.resize(nPartials);
	network.initPartialsIndices(reactionSize, reactionStartingIdx,
			reactionIndices);
	reactionVals.resize(nPartials);

	return;
}

void PetscSolver1DHandler::initializeConcentration(DM &da, Vec &C) {
	PetscErrorCode ierr;

	// Pointer for the concentration vector
	PetscScalar **concentrations = nullptr;
	ierr = DMDAVecGetArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolver1DHandler::initializeConcentration: "
			"DMDAVecGetArrayDOF failed.");

	// Get the local boundaries
	PetscInt xs, xm;
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr, "PetscSolver1DHandler::initializeConcentration: "
			"DMDAGetCorners failed.");

	// Get the last time step written in the HDF5 file
	int tempTimeStep = -2;
	bool hasConcentrations = false;
	if (!networkName.empty())
		hasConcentrations = xolotlCore::HDF5Utils::hasConcentrationGroup(
				networkName, tempTimeStep);

	// Give the surface position to the temperature handler
	temperatureHandler->updateSurfacePosition(grid[surfacePosition]);

	// Give the surface position to the temperature handler
	temperatureHandler->updateSurfacePosition(
			grid[surfacePosition + 1] - grid[1]);

	// Initialize the flux handler
	fluxHandler->initializeFluxHandler(network, surfacePosition, grid);

	// Initialize the grid for the diffusion
	diffusionHandler->initializeDiffusionGrid(advectionHandlers, grid);

	// Initialize the grid for the advection
	advectionHandlers[0]->initializeAdvectionGrid(advectionHandlers, grid);

	// Pointer for the concentration vector at a specific grid point
	PetscScalar *concOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	// + the super clusters
	const int dof = network.getDOF();

	// Get the single vacancy ID
	auto singleVacancyCluster = network.get(xolotlCore::Species::V, 1);
	int vacancyIndex = -1;
	if (singleVacancyCluster)
		vacancyIndex = singleVacancyCluster->getId() - 1;

	// Loop on all the grid points
	for (PetscInt i = xs; i < xs + xm; i++) {
		concOffset = concentrations[i];

		// Loop on all the clusters to initialize at 0.0
		for (int n = 0; n < dof - 1; n++) {
			concOffset[n] = 0.0;
		}

		// Temperature
		xolotlCore::Point<3> gridPosition { grid[i + 1] - grid[1], 0.0, 0.0 };
		concOffset[dof - 1] = temperatureHandler->getTemperature(gridPosition,
				0.0);

		// Initialize the vacancy concentration
		if (i >= surfacePosition + leftOffset && singleVacancyCluster
				&& !hasConcentrations && i < nX - rightOffset) {
			concOffset[vacancyIndex] = initialVConc;
		}
	}

	// If the concentration must be set from the HDF5 file
	if (hasConcentrations) {
		// Loop on the full grid
		for (int i = 0; i < nX; i++) {
			// Read the concentrations from the HDF5 file
			auto concVector = xolotlCore::HDF5Utils::readGridPoint(networkName,
					tempTimeStep, i);

			// Change the concentration only if we are on the locally owned part of the grid
			if (i >= xs && i < xs + xm) {
				concOffset = concentrations[i];
				// Loop on the concVector size
				for (unsigned int l = 0; l < concVector.size(); l++) {
					concOffset[(int) concVector.at(l).at(0)] =
							concVector.at(l).at(1);
				}
			}
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolver1DHandler::initializeConcentration: "
			"DMDAVecRestoreArrayDOF failed.");

	return;
}

void PetscSolver1DHandler::updateConcentration(TS &ts, Vec &localC, Vec &F,
		PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the local data vector from PETSc
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver1DHandler::updateConcentration: "
			"TSGetDM failed.");

	// Pointers to the PETSc arrays that start at the beginning (xs) of the
	// local array!
	PetscScalar **concs = nullptr, **updatedConcs = nullptr;
	// Get pointers to vector data
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver1DHandler::updateConcentration: "
			"DMDAVecGetArrayDOFRead (localC) failed.");
	ierr = DMDAVecGetArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolver1DHandler::updateConcentration: "
			"DMDAVecGetArrayDOF (F) failed.");

	// Get local grid boundaries
	PetscInt xs, xm;
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr, "PetscSolver1DHandler::updateConcentration: "
			"DMDAGetCorners failed.");

	// The following pointers are set to the first position in the conc or
	// updatedConc arrays that correspond to the beginning of the data for the
	// current grid point. They are accessed just like regular arrays.
	PetscScalar *concOffset = nullptr, *updatedConcOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network.getDOF();

	// Compute the total concentration of atoms contained in bubbles
	double atomConc = 0.0;

	// Loop over grid points to get the atom concentration
	// near the surface
	for (int xi = xs; xi < xs + xm; xi++) {
		// Boundary conditions
		if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset)
			continue;

		// We are only interested in the helium near the surface
		if (grid[xi] - grid[surfacePosition] > 2.0)
			continue;

		// Get the concentrations at this grid point
		concOffset = concs[xi];
		// Copy data into the PSIClusterReactionNetwork
		network.updateConcentrationsFromArray(concOffset);

		// Sum the total atom concentration
		atomConc += network.getTotalTrappedAtomConcentration()
				* (grid[xi + 1] - grid[xi]);
	}

	// Share the concentration with all the processes
	double totalAtomConc = 0.0;
	MPI_Allreduce(&atomConc, &totalAtomConc, 1, MPI_DOUBLE, MPI_SUM,
	MPI_COMM_WORLD);

	// Set the disappearing rate in the modified TM handler
	mutationHandler->updateDisappearingRate(totalAtomConc);

	// Declarations for variables used in the loop
	double **concVector = new double*[3];
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };

	// Loop over grid points computing ODE terms for each grid point
	for (PetscInt xi = xs; xi < xs + xm; xi++) {
		// Compute the old and new array offsets
		concOffset = concs[xi];
		updatedConcOffset = updatedConcs[xi];

		// Boundary conditions
		// Everything to the left of the surface is empty
		if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset) {
			continue;
		}

		// Set the grid position
		gridPosition[0] = grid[xi + 1] - grid[1];

		// Fill the concVector with the pointer to the middle, left, and right grid points
		concVector[0] = concOffset; // middle
		concVector[1] = concs[xi - 1]; // left
		concVector[2] = concs[xi + 1]; // right

		// Get the temperature from the temperature handler
		temperatureHandler->setTemperature(concOffset);
		double temperature = temperatureHandler->getTemperature(gridPosition,
				ftime);

		// Update the network if the temperature changed
		if (!xolotlCore::equal(temperature, lastTemperature)) {
			network.setTemperature(temperature);
			// Update the modified trap-mutation rate
			// that depends on the network reaction rates
			mutationHandler->updateTrapMutationRate(network);
			lastTemperature = temperature;
		}

		// Copy data into the ReactionNetwork so that it can
		// compute the fluxes properly. The network is only used to compute the
		// fluxes and hold the state data from the last time step. I'm reusing
		// it because it cuts down on memory significantly (about 400MB per
		// grid point) at the expense of being a little tricky to comprehend.
		network.updateConcentrationsFromArray(concOffset);

		// ----- Account for flux of incoming particles -----
		fluxHandler->computeIncidentFlux(ftime, updatedConcOffset, xi,
				surfacePosition);

		// ---- Compute the temperature over the locally owned part of the grid -----
		temperatureHandler->computeTemperature(concVector, updatedConcOffset,
				grid[xi + 1] - grid[xi], grid[xi + 2] - grid[xi + 1]);

		// ---- Compute diffusion over the locally owned part of the grid -----
		diffusionHandler->computeDiffusion(network, concVector,
				updatedConcOffset, grid[xi + 1] - grid[xi],
				grid[xi + 2] - grid[xi + 1], xi);

		// ---- Compute advection over the locally owned part of the grid -----
		for (int i = 0; i < advectionHandlers.size(); i++) {
			advectionHandlers[i]->computeAdvection(network, gridPosition,
					concVector, updatedConcOffset, grid[xi + 1] - grid[xi],
					grid[xi + 2] - grid[xi + 1], xi);
		}

		// ----- Compute the modified trap-mutation over the locally owned part of the grid -----
		mutationHandler->computeTrapMutation(network, concOffset,
				updatedConcOffset, xi);

		// ----- Compute the reaction fluxes over the locally owned part of the grid -----
		network.computeAllFluxes(updatedConcOffset);
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver1DHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOFRead (localC) failed.");
	ierr = DMDAVecRestoreArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolver1DHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOF (F) failed.");
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr, "PetscSolver1DHandler::updateConcentration: "
			"DMRestoreLocalVector failed.");

	// Clear memory
	delete[] concVector;

	return;
}

void PetscSolver1DHandler::computeOffDiagonalJacobian(TS &ts, Vec &localC,
		Mat &J, PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver1DHandler::computeOffDiagonalJacobian: "
			"TSGetDM failed.");

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network.getDOF();

	// Pointers to the PETSc arrays that start at the beginning (xs) of the
	// local array!
	PetscScalar **concs = nullptr;
	// Get pointers to vector data
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver1DHandler::computeOffDiagonalJacobian: "
			"DMDAVecGetArrayDOFRead (localC) failed.");

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset = nullptr;

	// Get local grid boundaries
	PetscInt xs, xm;
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr, "PetscSolver1DHandler::computeOffDiagonalJacobian: "
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
	MatStencil row, cols[3];
	PetscScalar diffVals[3 * nDiff];
	PetscInt diffIndices[nDiff];
	PetscScalar advecVals[2 * nAdvec];
	PetscInt advecIndices[nAdvec];
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };

	/*
	 Loop over grid points computing Jacobian terms for diffusion and advection
	 at each grid point
	 */
	for (PetscInt xi = xs; xi < xs + xm; xi++) {
		// Boundary conditions
		// Everything to the left of the surface is empty
		if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset)
			continue;

		// Set the grid position
		gridPosition[0] = grid[xi + 1] - grid[1];

		// Get the temperature from the temperature handler
		concOffset = concs[xi];
		temperatureHandler->setTemperature(concOffset);
		double temperature = temperatureHandler->getTemperature(gridPosition,
				ftime);

		// Update the network if the temperature changed
		if (!xolotlCore::equal(temperature, lastTemperature)) {
			network.setTemperature(temperature);
			// Update the modified trap-mutation rate
			// that depends on the network reaction rates
			mutationHandler->updateTrapMutationRate(network);
			lastTemperature = temperature;
		}

		// Get the partial derivatives for the temperature
		temperatureHandler->computePartialsForTemperature(diffVals, diffIndices,
				grid[xi + 1] - grid[xi], grid[xi + 2] - grid[xi + 1]);

		// Set grid coordinate and component number for the row
		row.i = xi;
		row.c = diffIndices[0];

		// Set grid coordinates and component numbers for the columns
		// corresponding to the middle, left, and right grid points
		cols[0].i = xi; // middle
		cols[0].c = diffIndices[0];
		cols[1].i = xi - 1; // left
		cols[1].c = diffIndices[0];
		cols[2].i = xi + 1; // right
		cols[2].c = diffIndices[0];

		ierr = MatSetValuesStencil(J, 1, &row, 3, cols, diffVals, ADD_VALUES);
		checkPetscError(ierr,
				"PetscSolver1DHandler::computeOffDiagonalJacobian: "
						"MatSetValuesStencil (temperature) failed.");

		// Get the partial derivatives for the diffusion
		diffusionHandler->computePartialsForDiffusion(network, diffVals,
				diffIndices, grid[xi + 1] - grid[xi],
				grid[xi + 2] - grid[xi + 1], xi);

		// Loop on the number of diffusion cluster to set the values in the Jacobian
		for (int i = 0; i < nDiff; i++) {
			// Set grid coordinate and component number for the row
			row.i = xi;
			row.c = diffIndices[i];

			// Set grid coordinates and component numbers for the columns
			// corresponding to the middle, left, and right grid points
			cols[0].i = xi; // middle
			cols[0].c = diffIndices[i];
			cols[1].i = xi - 1; // left
			cols[1].c = diffIndices[i];
			cols[2].i = xi + 1; // right
			cols[2].c = diffIndices[i];

			ierr = MatSetValuesStencil(J, 1, &row, 3, cols, diffVals + (3 * i),
					ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeOffDiagonalJacobian: "
							"MatSetValuesStencil (diffusion) failed.");
		}

		// Get the partial derivatives for the advection
		for (int l = 0; l < advectionHandlers.size(); l++) {
			advectionHandlers[l]->computePartialsForAdvection(network,
					advecVals, advecIndices, gridPosition,
					grid[xi + 1] - grid[xi], grid[xi + 2] - grid[xi + 1], xi);

			// Get the stencil indices to know where to put the partial derivatives in the Jacobian
			auto advecStencil = advectionHandlers[l]->getStencilForAdvection(
					gridPosition);

			// Get the number of advecting clusters
			nAdvec = advectionHandlers[l]->getNumberOfAdvecting();

			// Loop on the number of advecting cluster to set the values in the Jacobian
			for (int i = 0; i < nAdvec; i++) {
				// Set grid coordinate and component number for the row
				row.i = xi;
				row.c = advecIndices[i];

				// If we are on the sink, the partial derivatives are not the same
				// Both sides are giving their concentrations to the center
				if (advectionHandlers[l]->isPointOnSink(gridPosition)) {
					cols[0].i = xi - advecStencil[0]; // left?
					cols[0].c = advecIndices[i];
					cols[1].i = xi + advecStencil[0]; // right?
					cols[1].c = advecIndices[i];
				} else {
					// Set grid coordinates and component numbers for the columns
					// corresponding to the middle and other grid points
					cols[0].i = xi; // middle
					cols[0].c = advecIndices[i];
					cols[1].i = xi + advecStencil[0]; // left or right
					cols[1].c = advecIndices[i];
				}

				// Update the matrix
				ierr = MatSetValuesStencil(J, 1, &row, 2, cols,
						advecVals + (2 * i), ADD_VALUES);
				checkPetscError(ierr,
						"PetscSolver1DHandler::computeOffDiagonalJacobian: "
								"MatSetValuesStencil (advection) failed.");
			}
		}
	}

	return;
}

void PetscSolver1DHandler::computeDiagonalJacobian(TS &ts, Vec &localC, Mat &J,
		PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver1DHandler::computeDiagonalJacobian: "
			"TSGetDM failed.");

	// Get pointers to vector data
	PetscScalar **concs = nullptr;
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver1DHandler::computeDiagonalJacobian: "
			"DMDAVecGetArrayDOFRead failed.");

	// Get local grid boundaries
	PetscInt xs, xm;
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr, "PetscSolver1DHandler::computeDiagonalJacobian: "
			"DMDAGetCorners failed.");

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network.getDOF();

	// Compute the total concentration of atoms contained in bubbles
	double atomConc = 0.0;

	// Loop over grid points to get the atom concentration
	// near the surface
	for (int xi = xs; xi < xs + xm; xi++) {
		// Boundary conditions
		if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset)
			continue;

		// We are only interested in the helium near the surface
		if (grid[xi] - grid[surfacePosition] > 2.0)
			continue;

		// Get the concentrations at this grid point
		concOffset = concs[xi];
		// Copy data into the PSIClusterReactionNetwork
		network.updateConcentrationsFromArray(concOffset);

		// Sum the total atom concentration
		atomConc += network.getTotalTrappedAtomConcentration()
				* (grid[xi + 1] - grid[xi]);
	}

	// Share the concentration with all the processes
	double totalAtomConc = 0.0;
	MPI_Allreduce(&atomConc, &totalAtomConc, 1, MPI_DOUBLE, MPI_SUM,
	MPI_COMM_WORLD);

	// Set the disappearing rate in the modified TM handler
	mutationHandler->updateDisappearingRate(totalAtomConc);

	// Arguments for MatSetValuesStencil called below
	MatStencil rowId;
	MatStencil colIds[dof];
	int pdColIdsVectorSize = 0;

	// Declarations for variables used in the loop
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };

	// Loop over the grid points
	for (PetscInt xi = xs; xi < xs + xm; xi++) {
		// Boundary conditions
		// Everything to the left of the surface is empty
		if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset)
			continue;

		// Set the grid position
		gridPosition[0] = grid[xi + 1] - grid[1];

		// Get the temperature from the temperature handler
		concOffset = concs[xi];
		temperatureHandler->setTemperature(concOffset);
		double temperature = temperatureHandler->getTemperature(gridPosition,
				ftime);

		// Update the network if the temperature changed
		if (!xolotlCore::equal(temperature, lastTemperature)) {
			network.setTemperature(temperature);
			// Update the modified trap-mutation rate
			// that depends on the network reaction rates
			mutationHandler->updateTrapMutationRate(network);
			lastTemperature = temperature;
		}

		// Copy data into the ReactionNetwork so that it can
		// compute the new concentrations.
		network.updateConcentrationsFromArray(concOffset);

		// ----- Take care of the reactions for all the reactants -----

		// Compute all the partial derivatives for the reactions
		network.computeAllPartials(reactionStartingIdx, reactionIndices,
				reactionVals);

		// Update the column in the Jacobian that represents each DOF
		for (int i = 0; i < dof - 1; i++) {
			// Set grid coordinate and component number for the row
			rowId.i = xi;
			rowId.c = i;

			// Number of partial derivatives
			pdColIdsVectorSize = reactionSize[i];
			auto startingIdx = reactionStartingIdx[i];

			// Loop over the list of column ids
			for (int j = 0; j < pdColIdsVectorSize; j++) {
				// Set grid coordinate and component number for a column in the list
				colIds[j].i = xi;
				colIds[j].c = reactionIndices[startingIdx + j];
				// Get the partial derivative from the array of all of the partials
				reactingPartialsForCluster[j] = reactionVals[startingIdx + j];
			}
			// Update the matrix
			ierr = MatSetValuesStencil(J, 1, &rowId, pdColIdsVectorSize, colIds,
					reactingPartialsForCluster.data(), ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (reactions) failed.");
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
		int nMutating = mutationHandler->computePartialsForTrapMutation(network,
				mutationVals, mutationIndices, xi);

		// Loop on the number of helium undergoing trap-mutation to set the values
		// in the Jacobian
		for (int i = 0; i < nMutating; i++) {
			// Set grid coordinate and component number for the row and column
			// corresponding to the helium cluster
			row.i = xi;
			row.c = mutationIndices[3 * i];
			col.i = xi;
			col.c = mutationIndices[3 * i];

			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					mutationVals + (3 * i), ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (He trap-mutation) failed.");

			// Set component number for the row
			// corresponding to the HeV cluster created through trap-mutation
			row.c = mutationIndices[(3 * i) + 1];

			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					mutationVals + (3 * i) + 1, ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (HeV trap-mutation) failed.");

			// Set component number for the row
			// corresponding to the interstitial created through trap-mutation
			row.c = mutationIndices[(3 * i) + 2];

			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					mutationVals + (3 * i) + 2, ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (I trap-mutation) failed.");
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver1DHandler::computeDiagonalJacobian: "
			"DMDAVecRestoreArrayDOFRead failed.");
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr, "PetscSolver1DHandler::computeDiagonalJacobian: "
			"DMRestoreLocalVector failed.");

	return;
}

} /* end namespace xolotlSolver */
