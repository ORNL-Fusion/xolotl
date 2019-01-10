// Includes
#include <PetscSolver1DHandler.h>
#include <MathUtils.h>
#include <Constants.h>

namespace xcore = xolotlCore;

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

	// Get the MPI communicator on which to create the DMDA
	auto xolotlComm = xolotlCore::MPIUtils::getMPIComm();
	ierr = DMDACreate1d(xolotlComm, DM_BOUNDARY_MIRROR, nX, dof, 1,
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
	if (not networkName.empty() and movingSurface) {

		xolotlCore::XFile xfile(networkName);
		auto concGroup =
				xfile.getGroup<xolotlCore::XFile::ConcentrationGroup>();
		if (concGroup and concGroup->hasTimesteps()) {
			auto tsGroup = concGroup->getLastTimestepGroup();
			assert(tsGroup);
			surfacePosition = tsGroup->readSurface1D();
		}
	}

	// Initialize the surface of the first advection handler corresponding to the
	// advection toward the surface (or a dummy one if it is deactivated)
	advectionHandlers[0]->setLocation(grid[surfacePosition + 1] - grid[1]);

	// Prints the grid on one process
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);
	if (procId == 0) {
		for (int i = 1; i < grid.size() - 1; i++) {
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

	// Initialize the re-solution handler here
	// because it adds connectivity
	resolutionHandler->initialize(network);

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
	// Set it in the handler
	setLocalCoordinates(xs, xm);

	// Initialize the last temperature at each grid point on this process
	for (int i = 0; i < localXM + 2; i++) {
		lastTemperature.push_back(0.0);
	}
	network.addGridPoints(localXM + 2);

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
	temperatureHandler->updateSurfacePosition(surfacePosition);

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
	for (int i = localXS; i < localXS + localXM; i++) {
		concOffset = concentrations[i];

		// Loop on all the clusters to initialize at 0.0
		for (int n = 0; n < dof - 1; n++) {
			concOffset[n] = 0.0;
		}

		// Temperature
		xolotlCore::Point<3> gridPosition { (grid[i + 1]
				- grid[surfacePosition + 1])
				/ (grid[grid.size() - 2] - grid[surfacePosition + 1]), 0.0, 0.0 };
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
			// Set the temperature in the network
			double temp = myConcs[i][myConcs[i].size() - 1].second;
			network.setTemperature(temp, i);
			// Update the modified trap-mutation rate
			// that depends on the network reaction rates
			mutationHandler->updateTrapMutationRate(network);
			lastTemperature[i] = temp;
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolver1DHandler::initializeConcentration: "
			"DMDAVecRestoreArrayDOF failed.");

	// Set the rate for re-solution
	resolutionHandler->updateReSolutionRate(fluxHandler->getFluxAmplitude());

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
	for (int xi = localXS; xi < localXS + localXM; xi++) {
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
	auto xolotlComm = xolotlCore::MPIUtils::getMPIComm();
	MPI_Allreduce(&atomConc, &totalAtomConc, 1, MPI_DOUBLE, MPI_SUM,
			xolotlComm);

	// Set the disappearing rate in the modified TM handler
	mutationHandler->updateDisappearingRate(totalAtomConc);

	// Declarations for variables used in the loop
	double **concVector = new double*[3];
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };

	// Loop over grid points computing ODE terms for each grid point
	for (int xi = localXS; xi < localXS + localXM; xi++) {
		// Compute the old and new array offsets
		concOffset = concs[xi];
		updatedConcOffset = updatedConcs[xi];

		// Fill the concVector with the pointer to the middle, left, and right grid points
		concVector[0] = concOffset; // middle
		concVector[1] = concs[xi - 1]; // left
		concVector[2] = concs[xi + 1]; // right

		// Heat condition
		if (xi == surfacePosition) {
			temperatureHandler->computeTemperature(concVector,
					updatedConcOffset, grid[xi + 1] - grid[xi],
					grid[xi + 2] - grid[xi + 1], xi);
		}

		// Boundary conditions
		// Everything to the left of the surface is empty
		if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset) {
			continue;
		}
		// Free surface GB
		bool skip = false;
		for (auto &pair : gbVector) {
			if (xi == std::get<0>(pair)) {
				skip = true;
				break;
			}
		}
		if (skip)
			continue;

		// Update the network if the temperature changed
		// left
		double temperature = concs[xi - 1][dof - 1];
		if (std::fabs(lastTemperature[xi - localXS] - temperature) > 0.1) {
			network.setTemperature(temperature, xi - localXS);
			lastTemperature[xi - localXS] = temperature;
		}
		// right
		temperature = concs[xi + 1][dof - 1];
		if (std::fabs(lastTemperature[xi + 2 - localXS] - temperature) > 0.1) {
			network.setTemperature(temperature, xi + 2 - localXS);
			lastTemperature[xi + 2 - localXS] = temperature;
		}

		// Set the grid fraction
		gridPosition[0] = (grid[xi + 1] - grid[surfacePosition + 1])
				/ (grid[grid.size() - 1] - grid[surfacePosition + 1]);

		// Get the temperature from the temperature handler
		temperatureHandler->setTemperature(concOffset);
		temperature = temperatureHandler->getTemperature(gridPosition, ftime);
		// middle
		if (std::fabs(lastTemperature[xi + 1 - localXS] - temperature) > 0.1) {
			network.setTemperature(temperature, xi + 1 - localXS);
			// Update the modified trap-mutation rate
			// that depends on the network reaction rates
			mutationHandler->updateTrapMutationRate(network);
			lastTemperature[xi + 1 - localXS] = temperature;
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
				grid[xi + 1] - grid[xi], grid[xi + 2] - grid[xi + 1], xi);

		// ---- Compute diffusion over the locally owned part of the grid -----
		diffusionHandler->computeDiffusion(network, concVector,
				updatedConcOffset, grid[xi + 1] - grid[xi],
				grid[xi + 2] - grid[xi + 1], xi, localXS);

		// ---- Compute advection over the locally owned part of the grid -----
		// Set the grid position
		gridPosition[0] = grid[xi + 1] - grid[1];
		for (int i = 0; i < advectionHandlers.size(); i++) {
			advectionHandlers[i]->computeAdvection(network, gridPosition,
					concVector, updatedConcOffset, grid[xi + 1] - grid[xi],
					grid[xi + 2] - grid[xi + 1], xi, localXS);
		}

		// ----- Compute the modified trap-mutation over the locally owned part of the grid -----
		mutationHandler->computeTrapMutation(network, concOffset,
				updatedConcOffset, xi, localXS);

		// ----- Compute the re-solution over the locally owned part of the grid -----
		resolutionHandler->computeReSolution(network, concOffset,
				updatedConcOffset, xi, localXS);

		// ----- Compute the reaction fluxes over the locally owned part of the grid -----
		network.computeAllFluxes(updatedConcOffset, xi + 1 - localXS);
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

	// Pointers to the PETSc arrays that start at the beginning (xs) of the
	// local array!
	PetscScalar **concs = nullptr;
	// Get pointers to vector data
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver1DHandler::computeOffDiagonalJacobian: "
			"DMDAVecGetArrayDOFRead (localC) failed.");

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network.getDOF();

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
	for (int xi = localXS; xi < localXS + localXM; xi++) {
		// Heat condition
		if (xi == surfacePosition) {
			// Get the partial derivatives for the temperature
			temperatureHandler->computePartialsForTemperature(diffVals,
					diffIndices, grid[xi + 1] - grid[xi],
					grid[xi + 2] - grid[xi + 1], xi);

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

			ierr = MatSetValuesStencil(J, 1, &row, 3, cols, diffVals,
					ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeOffDiagonalJacobian: "
							"MatSetValuesStencil (temperature) failed.");
		}

		// Boundary conditions
		// Everything to the left of the surface is empty
		if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset)
			continue;
		// Free surface GB
		bool skip = false;
		for (auto &pair : gbVector) {
			if (xi == std::get<0>(pair)) {
				skip = true;
				break;
			}
		}
		if (skip)
			continue;

		// Update the network if the temperature changed
		// left
		double temperature = concs[xi - 1][dof - 1];
		if (std::fabs(lastTemperature[xi - localXS] - temperature) > 0.1) {
			network.setTemperature(temperature, xi - localXS);
			lastTemperature[xi - localXS] = temperature;
		}
		// right
		temperature = concs[xi + 1][dof - 1];
		if (std::fabs(lastTemperature[xi + 2 - localXS] - temperature) > 0.1) {
			network.setTemperature(temperature, xi + 2 - localXS);
			lastTemperature[xi + 2 - localXS] = temperature;
		}

		// Set the grid fraction
		gridPosition[0] = (grid[xi + 1] - grid[surfacePosition + 1])
				/ (grid[grid.size() - 1] - grid[surfacePosition + 1]);

		// Get the temperature from the temperature handler
		concOffset = concs[xi];
		temperatureHandler->setTemperature(concOffset);
		temperature = temperatureHandler->getTemperature(gridPosition, ftime);
		// middle
		if (std::fabs(lastTemperature[xi + 1 - localXS] - temperature) > 0.1) {
			network.setTemperature(temperature, xi + 1 - localXS);
			lastTemperature[xi + 1 - localXS] = temperature;
		}

		// Get the partial derivatives for the temperature
		temperatureHandler->computePartialsForTemperature(diffVals, diffIndices,
				grid[xi + 1] - grid[xi], grid[xi + 2] - grid[xi + 1], xi);

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
				grid[xi + 2] - grid[xi + 1], xi, localXS);

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
		// Set the grid position
		gridPosition[0] = grid[xi + 1] - grid[1];
		for (int l = 0; l < advectionHandlers.size(); l++) {
			advectionHandlers[l]->computePartialsForAdvection(network,
					advecVals, advecIndices, gridPosition,
					grid[xi + 1] - grid[xi], grid[xi + 2] - grid[xi + 1], xi,
					localXS);

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

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network.getDOF();

	// Compute the total concentration of atoms contained in bubbles
	double atomConc = 0.0;

	// Loop over grid points to get the atom concentration
	// near the surface
	for (int xi = localXS; xi < localXS + localXM; xi++) {
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
	auto xolotlComm = xolotlCore::MPIUtils::getMPIComm();
	MPI_Allreduce(&atomConc, &totalAtomConc, 1, MPI_DOUBLE, MPI_SUM,
			xolotlComm);

	// Set the disappearing rate in the modified TM handler
	mutationHandler->updateDisappearingRate(totalAtomConc);

	// Arguments for MatSetValuesStencil called below
	MatStencil rowId;
	MatStencil colIds[dof];
	int pdColIdsVectorSize = 0;

	// Declarations for variables used in the loop
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };

	// Loop over the grid points
	for (int xi = localXS; xi < localXS + localXM; xi++) {
		// Boundary conditions
		// Everything to the left of the surface is empty
		if (xi < surfacePosition + leftOffset || xi > nX - 1 - rightOffset)
			continue;
		// Free surface GB
		bool skip = false;
		for (auto &pair : gbVector) {
			if (xi == std::get<0>(pair)) {
				skip = true;
				break;
			}
		}
		if (skip)
			continue;

		// Set the grid fraction
		gridPosition[0] = (grid[xi + 1] - grid[surfacePosition + 1])
				/ (grid[grid.size() - 1] - grid[surfacePosition + 1]);

		// Get the temperature from the temperature handler
		concOffset = concs[xi];
		temperatureHandler->setTemperature(concOffset);
		double temperature = temperatureHandler->getTemperature(gridPosition,
				ftime);

		// Update the network if the temperature changed
		if (std::fabs(lastTemperature[xi + 1 - localXS] - temperature) > 0.1) {
			network.setTemperature(temperature, xi + 1 - localXS);
			// Update the modified trap-mutation rate
			// that depends on the network reaction rates
			mutationHandler->updateTrapMutationRate(network);
			lastTemperature[xi + 1 - localXS] = temperature;
		}

		// Copy data into the ReactionNetwork so that it can
		// compute the new concentrations.
		network.updateConcentrationsFromArray(concOffset);

		// ----- Take care of the reactions for all the reactants -----

		// Compute all the partial derivatives for the reactions
		network.computeAllPartials(reactionStartingIdx, reactionIndices,
				reactionVals, xi + 1 - localXS);

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
				mutationVals, mutationIndices, xi, localXS);

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

		// ----- Take care of the re-solution for all the reactants -----

		// Store the total number of Xe clusters in the network
		int nXenon = resolutionHandler->getNumberOfReSoluting();

		// Arguments for MatSetValuesStencil called below
		PetscScalar resolutionVals[10 * nXenon];
		PetscInt resolutionIndices[10 * nXenon];

		// Compute the partial derivative from re-solution at this grid point
		int nResoluting = resolutionHandler->computePartialsForReSolution(
				network, resolutionVals, resolutionIndices, xi, localXS);

		// Loop on the number of xenon to set the values in the Jacobian
		for (int i = 0; i < nResoluting; i++) {
			// Set grid coordinate and component number for the row and column
			// corresponding to the  large xenon cluster
			row.i = xi;
			row.c = resolutionIndices[10 * i];
			col.i = xi;
			col.c = resolutionIndices[10 * i];
			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					resolutionVals + (10 * i), ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (large Xe re-solution) failed.");
			col.c = resolutionIndices[(10 * i) + 1];
			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					resolutionVals + (10 * i) + 1, ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (large Xe re-solution) failed.");
			row.c = resolutionIndices[(10 * i) + 1];
			col.c = resolutionIndices[10 * i];
			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					resolutionVals + (10 * i) + 2, ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (large Xe re-solution) failed.");
			col.c = resolutionIndices[(10 * i) + 1];
			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					resolutionVals + (10 * i) + 3, ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (large Xe re-solution) failed.");

			// Set component number for the row
			// corresponding to the smaller xenon cluster created through re-solution
			row.c = resolutionIndices[(10 * i) + 4];
			col.c = resolutionIndices[10 * i];
			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					resolutionVals + (10 * i) + 4, ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (smaller Xe re-solution) failed.");
			col.c = resolutionIndices[(10 * i) + 1];
			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					resolutionVals + (10 * i) + 5, ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (smaller Xe re-solution) failed.");
			row.c = resolutionIndices[(10 * i) + 5];
			col.c = resolutionIndices[10 * i];
			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					resolutionVals + (10 * i) + 6, ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (smaller Xe re-solution) failed.");
			col.c = resolutionIndices[(10 * i) + 1];
			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					resolutionVals + (10 * i) + 7, ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (smaller Xe re-solution) failed.");

			// Set component number for the row
			// corresponding to the single xenon created through re-solution
			row.c = resolutionIndices[(10 * i) + 8];
			col.c = resolutionIndices[10 * i];
			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					resolutionVals + (10 * i) + 8, ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (Xe_1 re-solution) failed.");
			col.c = resolutionIndices[(10 * i) + 1];
			ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
					resolutionVals + (10 * i) + 9, ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolver1DHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (Xe_1 re-solution) failed.");
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
