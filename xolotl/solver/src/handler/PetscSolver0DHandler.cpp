// Includes
#include <xolotl/solver/handler/PetscSolver0DHandler.h>
#include <xolotl/util/MathUtils.h>
#include <xolotl/core/Constants.h>

namespace xolotl {
namespace solver {
namespace handler {

void PetscSolver0DHandler::createSolverContext(DM &da) {
	PetscErrorCode ierr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network.getDOF();

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, 1, dof + 1, 0,
	NULL, &da);
	checkPetscError(ierr, "PetscSolver0DHandler::createSolverContext: "
			"DMDACreate1d failed.");
	ierr = DMSetFromOptions(da);
	checkPetscError(ierr,
			"PetscSolver0DHandler::createSolverContext: DMSetFromOptions failed.");
	ierr = DMSetUp(da);
	checkPetscError(ierr,
			"PetscSolver0DHandler::createSolverContext: DMSetUp failed.");

	/*  The only spatial coupling in the Jacobian is due to diffusion.
	 *  The ofill (thought of as a dof by dof 2d (row-oriented) array represents
	 *  the nonzero coupling between degrees of freedom at one point with degrees
	 *  of freedom on the adjacent point to the left or right. A 1 at i,j in the
	 *  ofill array indicates that the degree of freedom i at a point is coupled
	 *  to degree of freedom j at the adjacent point.
	 *  In this case ofill has only a few diagonal entries since the only spatial
	 *  coupling is regular diffusion.
	 */
    core::network::IReactionNetwork::SparseFillMap ofill;

	// Initialize the temperature handler
	temperatureHandler->initializeTemperature(dof, ofill, dfill);

	// Initialize the nucleation handler here
	// because it adds connectivity
	nucleationHandler->initialize(network, dfill);

	// Get the diagonal fill
	auto nPartials = network.getDiagonalFill(dfill);

	// Load up the block fills
	auto dfillsparse = ConvertToPetscSparseFillMap(dof + 1, dfill);
	auto ofillsparse = ConvertToPetscSparseFillMap(dof + 1, ofill);
	ierr = DMDASetBlockFillsSparse(da, dfillsparse.data(), ofillsparse.data());
	checkPetscError(ierr, "PetscSolver0DHandler::createSolverContext: "
			"DMDASetBlockFills failed.");

	// Initialize the arrays for the reaction partial derivatives
	vals = Kokkos::View<double*>("solverPartials", nPartials);

	// Set the size of the partial derivatives vectors
	reactingPartialsForCluster.resize(dof, 0.0);

	return;
}

void PetscSolver0DHandler::initializeConcentration(DM &da, Vec &C) {
	PetscErrorCode ierr;

	// Initialize the last temperature
	temperature.push_back(0.0);

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
	auto singleVacancyCluster = network.getSingleVacancy();
	auto vacancyIndex = NetworkType::invalidIndex();
	if (singleVacancyCluster.getId() != NetworkType::invalidIndex())
		vacancyIndex = singleVacancyCluster.getId();

	// Get the concentration of the only grid point
	concOffset = concentrations[0];

	// Loop on all the clusters to initialize at 0.0
	for (int n = 0; n < dof; n++) {
		concOffset[n] = 0.0;
	}

	// Temperature
    plsm::SpaceVector<double, 3> gridPosition { 0.0, 0.0, 0.0 };
	concOffset[dof] = temperatureHandler->getTemperature(gridPosition, 0.0);
	temperature[0] = concOffset[dof];

	// Get the last time step written in the HDF5 file
	bool hasConcentrations = false;
	std::unique_ptr<io::XFile> xfile;
	std::unique_ptr<io::XFile::ConcentrationGroup> concGroup;
	if (not networkName.empty()) {

		xfile.reset(new io::XFile(networkName));
		concGroup = xfile->getGroup<io::XFile::ConcentrationGroup>();
		hasConcentrations = (concGroup and concGroup->hasTimesteps());
	}

	// Initialize the vacancy concentration
	if (vacancyIndex != NetworkType::invalidIndex() and not hasConcentrations) {
		concOffset[vacancyIndex] = initialVConc;
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

		for (auto const &currConcData : myConcs[0]) {
			concOffset[currConcData.first] = currConcData.second;
		}
		// Get the temperature
		double temp = myConcs[0][myConcs[0].size() - 1].second;
		temperature[0] = temp;
	}

	// Update the network with the temperature
	network.setTemperatures(temperature);
	network.syncClusterDataOnHost();

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolver0DHandler::initializeConcentration: "
			"DMDAVecRestoreArrayDOF failed.");

	// Set the rate for nucleation
	nucleationHandler->updateHeterogeneousNucleationRate(
			fluxHandler->getFluxAmplitude());

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
    plsm::SpaceVector<double, 3> gridPosition { 0.0, 0.0, 0.0 };

	// Get the old and new array offsets
	concOffset = concs[0];
	updatedConcOffset = updatedConcs[0];

	// Get the temperature from the temperature handler
	temperatureHandler->setTemperature(concOffset);
	double temp = temperatureHandler->getTemperature(gridPosition, ftime);

	// Update the network if the temperature changed
	if (std::fabs(temperature[0] - temp) > 0.1) {
		temperature[0] = temp;
		network.setTemperatures(temperature);
		network.syncClusterDataOnHost();
	}

	// ----- Account for flux of incoming particles -----
	fluxHandler->computeIncidentFlux(ftime, updatedConcOffset, 0, 0);

	// ----- Compute the heterogeneous nucleation -----
	nucleationHandler->computeHeterogeneousNucleation(network, concOffset,
			updatedConcOffset, 0, 0);

	// ----- Compute the reaction fluxes over the locally owned part of the grid -----
	using HostUnmanaged =
	Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concOffset, dof);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof);
	deep_copy(dConcs, hConcs);
	auto hFlux = HostUnmanaged(updatedConcOffset, dof);
	auto dFlux = Kokkos::View<double*>("Fluxes", dof);
	deep_copy(dFlux, hFlux);
	fluxCounter->increment();
	fluxTimer->start();
	network.computeAllFluxes(dConcs, dFlux, 0);
	fluxTimer->stop();
	deep_copy(hFlux, dFlux);

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver0DHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOFRead (localC) failed.");
	ierr = DMDAVecRestoreArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolver0DHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOF (F) failed.");

	return;
}

void PetscSolver0DHandler::computeJacobian(TS &ts, Vec &localC, Mat &J,
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
	MatStencil colId;
	int pdColIdsVectorSize = 0;

	// Set the grid position
    plsm::SpaceVector<double, 3> gridPosition { 0.0, 0.0, 0.0 };

	// Get the temperature from the temperature handler
	concOffset = concs[0];
	temperatureHandler->setTemperature(concOffset);
	double temp = temperatureHandler->getTemperature(gridPosition, ftime);

	// Update the network if the temperature changed
	if (std::fabs(temperature[0] - temp) > 0.1) {
		temperature[0] = temp;
		network.setTemperatures(temperature);
		network.syncClusterDataOnHost();
	}

	// ----- Take care of the reactions for all the reactants -----

	// Compute all the partial derivatives for the reactions
	using HostUnmanaged =
	Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concOffset, dof);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof);
	deep_copy(dConcs, hConcs);
	partialDerivativeCounter->increment();
	partialDerivativeTimer->start();
	network.computeAllPartials(dConcs, vals, 0);
	partialDerivativeTimer->stop();
	auto hPartials = create_mirror_view(vals);
	deep_copy(hPartials, vals);

	// Variable for the loop on reactants
	int startingIdx = 0;

	// Update the column in the Jacobian that represents each DOF
	for (int i = 0; i < dof; i++) {
		// Set grid coordinate and component number for the row
		rowId.i = 0;
		rowId.c = i;

		// Number of partial derivatives
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto &row = rowIter->second;
			pdColIdsVectorSize = row.size();

			// Loop over the list of column ids
			for (int j = 0; j < pdColIdsVectorSize; j++) {
				// Set grid coordinate and component number for a column in the list
				colIds[j].i = 0;
				colIds[j].c = row[j];
				// Get the partial derivative from the array of all of the partials
				reactingPartialsForCluster[j] = hPartials(startingIdx + j);
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

	// ----- Take care of the nucleation for all the reactants -----

	// Arguments for MatSetValuesStencil called below
	PetscScalar nucleationVals[2];
	PetscInt nucleationIndices[2];
	MatStencil rowIds[2];

	// Compute the partial derivative from nucleation at this grid point
	if (nucleationHandler->computePartialsForHeterogeneousNucleation(network,
			concOffset, nucleationVals, nucleationIndices, 0, 0)) {

		// Set grid coordinate and component number for the row and column
		// corresponding to the clusters involved in re-solution
		rowIds[0].i = 0;
		rowIds[0].c = nucleationIndices[0];
		rowIds[1].i = 0;
		rowIds[1].c = nucleationIndices[1];
		colIds[0].i = 0;
		colIds[0].c = nucleationIndices[0];
		ierr = MatSetValuesStencil(J, 2, rowIds, 1, colIds, nucleationVals,
				ADD_VALUES);
		checkPetscError(ierr, "PetscSolver0DHandler::computeDiagonalJacobian: "
				"MatSetValuesStencil (Xe nucleation) failed.");
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver0DHandler::computeDiagonalJacobian: "
			"DMDAVecRestoreArrayDOFRead failed.");

	return;
}

} /* end namespace handler */
} /* end namespace solver */
} /* end namespace xolotl */
