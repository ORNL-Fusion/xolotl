// Includes
#include <PetscSolverExpHandler.h>
#include <MathUtils.h>
#include <Constants.h>

namespace xolotlSolver {

template<typename TImpl>
void PetscSolverExpHandler<TImpl>::createSolverContext(DM &da) {
	PetscErrorCode ierr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = expNetwork.getDOF();

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, 1, dof, 0,
	NULL, &da);
	checkPetscError(ierr, "PetscSolverExpHandler::createSolverContext: "
			"DMDACreate1d failed.");
	ierr = DMSetFromOptions(da);
	checkPetscError(ierr,
			"PetscSolverExpHandler::createSolverContext: DMSetFromOptions failed.");
	ierr = DMSetUp(da);
	checkPetscError(ierr,
			"PetscSolverExpHandler::createSolverContext: DMSetUp failed.");

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
	temperatureHandler->initializeTemperature(network, ofill, dfill);

	// Get the diagonal fill
	auto nPartials = expNetwork.getDiagonalFill(dfill);

	// Load up the block fills
	auto dfillsparse = ConvertToPetscSparseFillMap(dof, dfill);
	auto ofillsparse = ConvertToPetscSparseFillMap(dof, ofill);
	ierr = DMDASetBlockFillsSparse(da, dfillsparse.data(), ofillsparse.data());
	checkPetscError(ierr, "PetscSolverExpHandler::createSolverContext: "
			"DMDASetBlockFills failed.");

	// Initialize the arrays for the reaction partial derivatives
	expVals = Kokkos::View<double*, Kokkos::MemoryUnmanaged>(nPartials);

	// Set the size of the partial derivatives vectors
	reactingPartialsForCluster.resize(dof, 0.0);

	return;
}

template<typename TImpl>
void PetscSolverExpHandler<TImpl>::initializeConcentration(DM &da, Vec &C) {
	PetscErrorCode ierr;

	// Initialize the last temperature and rates
	lastTemperature.push_back(0.0);
	network.addGridPoints(1);

	// Pointer for the concentration vector
	PetscScalar **concentrations = nullptr;
	ierr = DMDAVecGetArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolverExpHandler::initializeConcentration: "
			"DMDAVecGetArrayDOF failed.");

	// Initialize the flux handler
	fluxHandler->initializeFluxHandler(network, 0, grid);

	// Pointer for the concentration vector at a specific grid point
	PetscScalar *concOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	// + the super clusters
	const int dof = expNetwork.getDOF();

	// Get the concentration of the only grid point
	concOffset = concentrations[0];

	// Loop on all the clusters to initialize at 0.0
	for (int n = 0; n < dof - 1; n++) {
		concOffset[n] = 0.0;
	}

	// Temperature
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };
	concOffset[dof - 1] = temperatureHandler->getTemperature(gridPosition, 0.0);

	// Get the last time step written in the HDF5 file
	bool hasConcentrations = false;
	std::unique_ptr<xolotlCore::XFile> xfile;
	std::unique_ptr<xolotlCore::XFile::ConcentrationGroup> concGroup;
	if (not networkName.empty()) {

		xfile.reset(new xolotlCore::XFile(networkName));
		concGroup = xfile->getGroup<xolotlCore::XFile::ConcentrationGroup>();
		hasConcentrations = (concGroup and concGroup->hasTimesteps());
	}

	// If the concentration must be set from the HDF5 file
	if (hasConcentrations) {
		// Read the concentrations from the HDF5 file for
		// each of our grid points.
		assert(concGroup);
		auto tsGroup = concGroup->getLastTimestepGroup();
		assert(tsGroup);
		auto myConcs = tsGroup->readConcentrations(*xfile, 0, 1);

		// Apply the concentrations we just read.
		concOffset = concentrations[0];

		for (auto const& currConcData : myConcs[0]) {
			concOffset[currConcData.first] = currConcData.second;
		}
		// Set the temperature in the network
		double temp = myConcs[0][myConcs[0].size() - 1].second;
		network.setTemperature(temp, 0);
		lastTemperature[0] = temp;
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolverExpHandler::initializeConcentration: "
			"DMDAVecRestoreArrayDOF failed.");

	return;
}

template<typename TImpl>
void PetscSolverExpHandler<TImpl>::updateConcentration(TS &ts, Vec &localC,
		Vec &F, PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the local data vector from PETSc
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolverExpHandler::updateConcentration: "
			"TSGetDM failed.");

	// Pointers to the PETSc arrays that start at the beginning (xs) of the
	// local array!
	PetscScalar **concs = nullptr, **updatedConcs = nullptr;
	// Get pointers to vector data
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolverExpHandler::updateConcentration: "
			"DMDAVecGetArrayDOFRead (localC) failed.");
	ierr = DMDAVecGetArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolverExpHandler::updateConcentration: "
			"DMDAVecGetArrayDOF (F) failed.");

	// The following pointers are set to the first position in the conc or
	// updatedConc arrays that correspond to the beginning of the data for the
	// current grid point. They are accessed just like regular arrays.
	PetscScalar *concOffset = nullptr, *updatedConcOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = expNetwork.getDOF();

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
	if (std::fabs(lastTemperature[0] - temperature) > 0.1) {
		network.setTemperature(temperature);
		lastTemperature[0] = temperature;
	}

	// ----- Account for flux of incoming particles -----
	fluxHandler->computeIncidentFlux(ftime, updatedConcOffset, 0, 0);

	// ----- Compute the reaction fluxes over the locally owned part of the grid -----
	ConcentrationsView xView(concOffset, dof);
	FluxesView fView(updatedConcOffset, dof);
	expNetwork.computeAllFluxes(xView, fView, 0);

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolverExpHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOFRead (localC) failed.");
	ierr = DMDAVecRestoreArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolverExpHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOF (F) failed.");

	return;
}

template<typename TImpl>
void PetscSolverExpHandler<TImpl>::computeOffDiagonalJacobian(TS &ts,
		Vec &localC, Mat &J, PetscReal ftime) {
	// Does nothing in 0D

	return;
}

template<typename TImpl>
void PetscSolverExpHandler<TImpl>::computeDiagonalJacobian(TS &ts, Vec &localC,
		Mat &J, PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolverExpHandler::computeDiagonalJacobian: "
			"TSGetDM failed.");

	// Get pointers to vector data
	PetscScalar **concs = nullptr;
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolverExpHandler::computeDiagonalJacobian: "
			"DMDAVecGetArrayDOFRead failed.");

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = expNetwork.getDOF();

	// Arguments for MatSetValuesStencil called below
	MatStencil rowId;
	MatStencil colIds[dof];
	MatStencil colId;
	int pdColIdsVectorSize = 0;

	// Set the grid position
	xolotlCore::Point<3> gridPosition { 0.0, 0.0, 0.0 };

	// Get the temperature from the temperature handler
	concOffset = concs[0];
	temperatureHandler->setTemperature(concOffset);
	double temperature = temperatureHandler->getTemperature(gridPosition,
			ftime);

	// Update the network if the temperature changed
	if (std::fabs(lastTemperature[0] - temperature) > 0.1) {
		network.setTemperature(temperature);
		lastTemperature[0] = temperature;
	}

	// ----- Take care of the reactions for all the reactants -----

	// Compute all the partial derivatives for the reactions
	ConcentrationsView xView(concOffset, dof);
	expNetwork.computeAllPartials(concOffset, expVals, 0);

	// Variable for the loop on reactants
	int startingIdx = 0;
	// Update the column in the Jacobian that represents each DOF
	for (int i = 0; i < dof - 1; i++) {
		// Set grid coordinate and component number for the row
		rowId.i = 0;
		rowId.c = i;

		// Number of partial derivatives
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
            const auto& row = rowIter->second;
			pdColIdsVectorSize = row.size();

			// Loop over the list of column ids
			for (int j = 0; j < pdColIdsVectorSize; j++) {
				// Set grid coordinate and component number for a column in the list
				colIds[j].i = 0;
				colIds[j].c = row[j];
				// Get the partial derivative from the array of all of the partials
				reactingPartialsForCluster[j] = expVals(startingIdx + j);
			}
			// Update the matrix
			ierr = MatSetValuesStencil(J, 1, &rowId, pdColIdsVectorSize, colIds,
					reactingPartialsForCluster.data(), ADD_VALUES);
			checkPetscError(ierr,
					"PetscSolverExpHandler::computeDiagonalJacobian: "
							"MatSetValuesStencil (reactions) failed.");

			// Increase the starting index
			startingIdx += pdColIdsVectorSize;
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolverExpHandler::computeDiagonalJacobian: "
			"DMDAVecRestoreArrayDOFRead failed.");

	return;
}

} /* end namespace xolotlSolver */
