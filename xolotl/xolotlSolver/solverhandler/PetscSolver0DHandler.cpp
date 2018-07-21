// Includes
#include <PetscSolver0DHandler.h>
#include <MathUtils.h>
#include <Constants.h>

namespace xolotlSolver {

void PetscSolver0DHandler::createSolverContext(DM &da) {
	PetscErrorCode ierr;
	// Recompute Ids and network size and redefine the connectivities
	network.reinitializeConnectivities();

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network.getDOF();

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, 1, dof, 0,
	NULL, &da);
	checkPetscError(ierr, "PetscSolver0DHandler::createSolverContext: "
			"DMDACreate1d failed.");
	ierr = DMSetFromOptions(da);
	checkPetscError(ierr,
			"PetscSolver0DHandler::createSolverContext: DMSetFromOptions failed.");
	ierr = DMSetUp(da);
	checkPetscError(ierr,
			"PetscSolver0DHandler::createSolverContext: DMSetUp failed.");

	// Set the size of the partial derivatives vectors
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

	// Get the diagonal fill
	network.getDiagonalFill(dfill);

	// Load up the block fills
	auto dfillsparse = ConvertToPetscSparseFillMap(dof, dfill);
	auto ofillsparse = ConvertToPetscSparseFillMap(dof, ofill);
	ierr = DMDASetBlockFillsSparse(da, dfillsparse.data(), ofillsparse.data());
	checkPetscError(ierr, "PetscSolver0DHandler::createSolverContext: "
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

void PetscSolver0DHandler::initializeConcentration(DM &da, Vec &C) {
	PetscErrorCode ierr;

	// Pointer for the concentration vector
	PetscScalar **concentrations = nullptr;
	ierr = DMDAVecGetArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolver0DHandler::initializeConcentration: "
			"DMDAVecGetArrayDOF failed.");

	// Initialize the flux handler
	fluxHandler->initializeFluxHandler(network, 0, grid);

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

	// Get the concentration of the only grid point
	concOffset = concentrations[0];

	// Loop on all the clusters to initialize at 0.0
	for (int n = 0; n < dof - 1; n++) {
		concOffset[n] = 0.0;
	}

	// Temperature
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };
	concOffset[dof - 1] = temperatureHandler->getTemperature(gridPosition, 0.0);

    // Determine if the HDF5 file has any concentrations.
    // TODO are we assuming networkName is non-empty here?
    xolotlCore::XFile xfile(networkName);
    auto concGroup = xfile.getGroup<xolotlCore::XFile::ConcentrationGroup>();
    bool hasConcentrations = (concGroup and concGroup->hasTimesteps());

	// Initialize the vacancy concentration
	if (singleVacancyCluster and not hasConcentrations) {
		concOffset[vacancyIndex] = initialVConc;
	}

	// If the concentration must be set from the HDF5 file
	if (hasConcentrations) {
		// Read the concentrations from the HDF5 file
        auto tsGroup = concGroup->getLastTimestepGroup();
        auto concVector = tsGroup->readGridPoint(0);

		concOffset = concentrations[0];
		// Loop on the concVector size
		for (unsigned int l = 0; l < concVector.size(); l++) {
			concOffset[(int) concVector.at(l).at(0)] = concVector.at(l).at(1);
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolver0DHandler::initializeConcentration: "
			"DMDAVecRestoreArrayDOF failed.");

	return;
}

void PetscSolver0DHandler::updateConcentration(TS &ts, Vec &localC, Vec &F,
		PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the local data vector from PETSc
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver0DHandler::updateConcentration: "
			"TSGetDM failed.");

	// Pointers to the PETSc arrays that start at the beginning (xs) of the
	// local array!
	PetscScalar **concs = nullptr, **updatedConcs = nullptr;
	// Get pointers to vector data
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver0DHandler::updateConcentration: "
			"DMDAVecGetArrayDOFRead (localC) failed.");
	ierr = DMDAVecGetArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolver0DHandler::updateConcentration: "
			"DMDAVecGetArrayDOF (F) failed.");

	// The following pointers are set to the first position in the conc or
	// updatedConc arrays that correspond to the beginning of the data for the
	// current grid point. They are accessed just like regular arrays.
	PetscScalar *concOffset = nullptr, *updatedConcOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network.getDOF();

	// Set the grid position
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };

	// Get the old and new array offsets
	concOffset = concs[0];
	updatedConcOffset = updatedConcs[0];

	// Get the temperature from the temperature handler
	temperatureHandler->setTemperature(concOffset);
	double temperature = temperatureHandler->getTemperature(gridPosition,
			ftime);

	// Update the network if the temperature changed
	if (!xolotlCore::equal(temperature, lastTemperature)) {
		network.setTemperature(temperature);
		lastTemperature = temperature;
	}

	// Copy data into the ReactionNetwork so that it can
	// compute the fluxes properly. The network is only used to compute the
	// fluxes and hold the state data from the last time step. I'm reusing
	// it because it cuts down on memory significantly (about 400MB per
	// grid point) at the expense of being a little tricky to comprehend.
	network.updateConcentrationsFromArray(concOffset);

	// ----- Account for flux of incoming particles -----
	fluxHandler->computeIncidentFlux(ftime, updatedConcOffset, 0, 0);

	// ----- Compute the reaction fluxes over the locally owned part of the grid -----
	network.computeAllFluxes(updatedConcOffset);

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver0DHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOFRead (localC) failed.");
	ierr = DMDAVecRestoreArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolver0DHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOF (F) failed.");
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr, "PetscSolver0DHandler::updateConcentration: "
			"DMRestoreLocalVector failed.");

	return;
}

void PetscSolver0DHandler::computeOffDiagonalJacobian(TS &ts, Vec &localC,
		Mat &J, PetscReal ftime) {
	// Does nothing in 0D

	return;
}

void PetscSolver0DHandler::computeDiagonalJacobian(TS &ts, Vec &localC, Mat &J,
		PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver0DHandler::computeDiagonalJacobian: "
			"TSGetDM failed.");

	// Get pointers to vector data
	PetscScalar **concs = nullptr;
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver0DHandler::computeDiagonalJacobian: "
			"DMDAVecGetArrayDOFRead failed.");

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network.getDOF();

	// Arguments for MatSetValuesStencil called below
	MatStencil rowId;
	MatStencil colIds[dof];
	int pdColIdsVectorSize = 0;

	// Set the grid position
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };

	// Get the temperature from the temperature handler
	concOffset = concs[0];
	temperatureHandler->setTemperature(concOffset);
	double temperature = temperatureHandler->getTemperature(gridPosition,
			ftime);

	// Update the network if the temperature changed
	if (!xolotlCore::equal(temperature, lastTemperature)) {
		network.setTemperature(temperature);
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
		rowId.i = 0;
		rowId.c = i;

		// Number of partial derivatives
		pdColIdsVectorSize = reactionSize[i];
		auto startingIdx = reactionStartingIdx[i];

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set grid coordinate and component number for a column in the list
			colIds[j].i = 0;
			colIds[j].c = reactionIndices[startingIdx + j];
			// Get the partial derivative from the array of all of the partials
			reactingPartialsForCluster[j] = reactionVals[startingIdx + j];
		}
		// Update the matrix
		ierr = MatSetValuesStencil(J, 1, &rowId, pdColIdsVectorSize, colIds,
				reactingPartialsForCluster.data(), ADD_VALUES);
		checkPetscError(ierr, "PetscSolver0DHandler::computeDiagonalJacobian: "
				"MatSetValuesStencil (reactions) failed.");
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver0DHandler::computeDiagonalJacobian: "
			"DMDAVecRestoreArrayDOFRead failed.");
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr, "PetscSolver0DHandler::computeDiagonalJacobian: "
			"DMRestoreLocalVector failed.");

	return;
}

} /* end namespace xolotlSolver */
