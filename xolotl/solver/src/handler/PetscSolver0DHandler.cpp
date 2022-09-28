#include <petscconf.h>
#include <petscdmda_kokkos.hpp>

#include <xolotl/core/Constants.h>
#include <xolotl/core/Types.h>
#include <xolotl/io/XFile.h>
#include <xolotl/solver/handler/PetscSolver0DHandler.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace solver
{
namespace handler
{
void
PetscSolver0DHandler::createSolverContext(DM& da)
{
	PetscErrorCode ierr;

	// Degrees of freedom is the total number of clusters in the network
	// + moments
	const auto dof = network.getDOF();

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	XOLOTL_LOG << "SolverHandler: 0D simulation";
	for (auto pair : initialConc) {
		XOLOTL_LOG << ", initial concentration for Id: " << pair.first
				   << " of: " << pair.second << " nm-3";
	}

	// Get the MPI communicator on which to create the DMDA
	auto xolotlComm = util::getMPIComm();
	ierr = DMDACreate1d(xolotlComm, DM_BOUNDARY_NONE, 1, dof + 1, 0, NULL, &da);
	checkPetscError(ierr,
		"PetscSolver0DHandler::createSolverContext: "
		"DMDACreate1d failed.");
	ierr = DMSetFromOptions(da);
	checkPetscError(ierr,
		"PetscSolver0DHandler::createSolverContext: DMSetFromOptions failed.");
	ierr = DMSetUp(da);
	checkPetscError(
		ierr, "PetscSolver0DHandler::createSolverContext: DMSetUp failed.");
}

void
PetscSolver0DHandler::initializeSolverContext(DM& da, TS& ts)
{
	PetscErrorCode ierr;

	// Degrees of freedom is the total number of clusters in the network
	// + moments
	const auto dof = network.getDOF();

	/* The ofill (thought of as a dof by dof 2d (row-oriented) array represents
	 * the nonzero coupling between degrees of freedom at one point with
	 * degrees of freedom on the adjacent point to the left or right.
	 */
	core::network::IReactionNetwork::SparseFillMap ofill;

	// Initialize the temperature handler
	temperatureHandler->initialize(dof);

	// Get the diagonal fill
	auto nPartials = network.getDiagonalFill(dfill);

	// Preallocate matrix
	Mat J;
	ierr = TSGetRHSJacobian(ts, &J, nullptr, nullptr, nullptr);
	checkPetscError(ierr,
		"PetscSolver0DHandler::initializeSolverContext: "
		"TSGetRHSJacobian failed.");
	auto [rows, cols] = convertToCoordinateListPair(dof, dfill);
    // handling temperature (FIXME)
	rows.push_back(dof);
	cols.push_back(dof);
    ++nPartials;
    //
	ierr = MatSetPreallocationCOO(J, rows.size(), rows.data(), cols.data());
	checkPetscError(ierr,
		"PetscSolver0DHandler::initializeSolverContext: "
		"MatSetPreallocationCOO failed.");

	// Initialize the arrays for the reaction partial derivatives
	vals = Kokkos::View<double*>("solverPartials", nPartials);

	// Set the size of the partial derivatives vectors
	reactingPartialsForCluster.resize(dof, 0.0);

	// Initialize the flux handler
	fluxHandler->initializeFluxHandler(network, 0, grid);
}

void
PetscSolver0DHandler::initializeConcentration(DM& da, Vec& C)
{
	PetscErrorCode ierr;

	// Initialize the last temperature
	temperature.push_back(0.0);

	// Pointer for the concentration vector
	PetscScalar** concentrations = nullptr;
	ierr = DMDAVecGetArrayDOF(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver0DHandler::initializeConcentration: "
		"DMDAVecGetArrayDOF failed.");

	// Pointer for the concentration vector at a specific grid point
	PetscScalar* concOffset = nullptr;

	// Degrees of freedom is the total number of clusters in the network
	// + moments
	const auto dof = network.getDOF();

	// Get the concentration of the only grid point
	concOffset = concentrations[0];

	// Loop on all the clusters to initialize at 0.0
	for (auto n = 0; n < dof; n++) {
		concOffset[n] = 0.0;
	}

	// Temperature
	plsm::SpaceVector<double, 3> gridPosition{0.0, 0.0, 0.0};
	concOffset[dof] = temperatureHandler->getTemperature(gridPosition, 0.0);
	temperature[0] = concOffset[dof];

	// Get the last time step written in the HDF5 file
	bool hasConcentrations = false;
	std::unique_ptr<io::XFile> xfile;
	std::unique_ptr<io::XFile::ConcentrationGroup> concGroup;
	if (not networkName.empty()) {
		xfile = std::make_unique<io::XFile>(networkName);
		concGroup = xfile->getGroup<io::XFile::ConcentrationGroup>();
		hasConcentrations = (concGroup and concGroup->hasTimesteps());
	}

	// Initialize the option specified concentration
	if (hasConcentrations) {
		for (auto pair : initialConc) {
			concOffset[pair.first] = pair.second;
		}
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
		// Get the temperature
		double temp = myConcs[0][myConcs[0].size() - 1].second;
		temperature[0] = temp;
	}

	// Update the network with the temperature
	auto depths = std::vector<double>(1, 1.0);
	network.setTemperatures(temperature, depths);

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver0DHandler::initializeConcentration: "
		"DMDAVecRestoreArrayDOF failed.");

	return;
}

std::vector<std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
PetscSolver0DHandler::getConcVector(DM& da, Vec& C)
{
	// Initial declaration
	PetscErrorCode ierr;
	const double* gridPointSolution = nullptr;

	// Pointer for the concentration vector
	PetscScalar** concentrations = nullptr;
	ierr = DMDAVecGetArrayDOFRead(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver0DHandler::getConcVector: "
		"DMDAVecGetArrayDOFRead failed.");

	// Get the network and dof
	auto& network = getNetwork();
	const auto dof = network.getDOF();

	// Create the vector for the concentrations
	std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
		toReturn;

	// Access the solution data for the current grid point.
	gridPointSolution = concentrations[0];

	// Create the temporary vector for this grid point
	std::vector<std::pair<IdType, double>> tempVector;
	for (auto l = 0; l < dof + 1; ++l) {
		if (std::fabs(gridPointSolution[l]) > 1.0e-16) {
			tempVector.push_back(std::make_pair(l, gridPointSolution[l]));
		}
	}
	std::vector<std::vector<std::pair<IdType, double>>> tempTempVector;
	tempTempVector.push_back(tempVector);
	std::vector<std::vector<std::vector<std::pair<IdType, double>>>>
		tempTempTempVector;
	tempTempTempVector.push_back(tempTempVector);
	toReturn.push_back(tempTempTempVector);

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver0DHandler::getConcVector: "
		"DMDAVecRestoreArrayDOFRead failed.");

	return toReturn;
}

void
PetscSolver0DHandler::setConcVector(DM& da, Vec& C,
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
		"PetscSolver0DHandler::setConcVector: "
		"DMDAVecGetArrayDOF failed.");

	// Get the DOF of the network
	const auto dof = network.getDOF();

	// Get the local concentration
	gridPointSolution = concentrations[0];

	// Loop on the given vector
	for (auto l = 0; l < concVector[0][0][0].size(); l++) {
		gridPointSolution[concVector[0][0][0][l].first] =
			concVector[0][0][0][l].second;
	}

	// Set the temperature in the network
	temperature[0] = gridPointSolution[dof];
	auto depths = std::vector<double>(1, 1.0);
	network.setTemperatures(temperature, depths);

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr,
		"PetscSolver0DHandler::setConcVector: "
		"DMDAVecRestoreArrayDOF failed.");

	return;
}

void
PetscSolver0DHandler::updateConcentration(
	TS& ts, Vec& localC, Vec& F, PetscReal ftime)
{
	PetscErrorCode ierr;

	// Get the local data vector from PETSc
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr,
		"PetscSolver0DHandler::updateConcentration: "
		"TSGetDM failed.");

	// Pointers to the PETSc arrays that start at the beginning of the
	// local array
	PetscOffsetView<const PetscScalar**> concs;
	ierr = DMDAVecGetKokkosOffsetViewDOF(da, localC, &concs);
	checkPetscError(ierr,
		"PetscSolver0DHandler::updateConcentration: "
		"DMDAVecGetKokkosOffsetViewDOF (localC) failed.");
	PetscOffsetView<PetscScalar**> updatedConcs;
	ierr = DMDAVecGetKokkosOffsetViewDOFWrite(da, F, &updatedConcs);
	checkPetscError(ierr,
		"PetscSolver0DHandler::updateConcentration: "
		"DMDAVecGetKokkosOffsetViewDOFWrite (F) failed.");

	// Set the grid position
	plsm::SpaceVector<double, 3> gridPosition{0.0, 0.0, 0.0};

	// The following pointers are set to the first position in the conc or
	// updatedConc arrays that correspond to the beginning of the data for the
	// current grid point.
	auto concOffset = subview(concs, 0, Kokkos::ALL).view();
	auto updatedConcOffset = subview(updatedConcs, 0, Kokkos::ALL).view();

	// Degrees of freedom is the total number of clusters in the network +
	// moments
	const auto dof = network.getDOF();

	// Update the time in the network
	network.setTime(ftime);

	// Get the temperature from the temperature handler
	temperatureHandler->setTemperature(concOffset);
	double temp = temperatureHandler->getTemperature(gridPosition, ftime);

	// Update the network if the temperature changed
	if (std::fabs(temperature[0] - temp) > 0.1) {
		temperature[0] = temp;
		auto depths = std::vector<double>(1, 1.0);
		network.setTemperatures(temperature, depths);
	}

	// ----- Account for flux of incoming particles -----
	fluxHandler->computeIncidentFlux(ftime, updatedConcOffset, 0, 0);

	// ----- Compute the reaction fluxes over the locally owned part of the grid
	// -----
	fluxCounter->increment();
	fluxTimer->start();
	network.computeAllFluxes(concOffset, updatedConcOffset);
	fluxTimer->stop();

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreKokkosOffsetViewDOF(da, localC, &concs);
	checkPetscError(ierr,
		"PetscSolver0DHandler::updateConcentration: "
		"DMDAVecRestoreKokkosOffsetViewDOF (localC) failed.");
	ierr = DMDAVecRestoreKokkosOffsetViewDOFWrite(da, F, &updatedConcs);
	checkPetscError(ierr,
		"PetscSolver0DHandler::updateConcentration: "
		"DMDAVecRestoreKokkosOffsetViewDOFWrite (F) failed.");
}

void
PetscSolver0DHandler::computeJacobian(
	TS& ts, Vec& localC, Mat& J, PetscReal ftime)
{
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr,
		"PetscSolver0DHandler::computeDiagonalJacobian: "
		"TSGetDM failed.");

	// Get pointers to vector data
	PetscOffsetView<const PetscScalar**> concs;
	ierr = DMDAVecGetKokkosOffsetViewDOF(da, localC, &concs);
	checkPetscError(ierr,
		"PetscSolver0DHandler::computeDiagonalJacobian: "
		"DMDAVecGetKokkosOffsetViewDOF failed.");

	// Degrees of freedom is the total number of clusters in the network +
	// moments
	const auto dof = network.getDOF();

	// Arguments for MatSetValuesStencil called below
	MatStencil rowId;
	MatStencil colIds[dof];
	MatStencil colId;
	IdType pdColIdsVectorSize = 0;

	// Set the grid position
	plsm::SpaceVector<double, 3> gridPosition{0.0, 0.0, 0.0};

	// Update the time in the network
	network.setTime(ftime);

	// Get the temperature from the temperature handler
	auto concOffset = subview(concs, 0, Kokkos::ALL).view();
	temperatureHandler->setTemperature(concOffset);
	double temp = temperatureHandler->getTemperature(gridPosition, ftime);

	// Update the network if the temperature changed
	if (std::fabs(temperature[0] - temp) > 0.1) {
		temperature[0] = temp;
		auto depths = std::vector<double>(1, 1.0);
		network.setTemperatures(temperature, depths);
	}

	// ----- Take care of the reactions for all the reactants -----

	// Compute all the partial derivatives for the reactions
	partialDerivativeCounter->increment();
	partialDerivativeTimer->start();
	network.computeAllPartials(concOffset, vals);
	partialDerivativeTimer->stop();

	ierr = MatSetValuesCOO(J, vals.data(), ADD_VALUES);
	checkPetscError(
		ierr, "PetscSolver0DHandler::computeJacobian: MatSetValuesCOO failed.");

	// Reset the values
    resetJacobianValues();

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreKokkosOffsetViewDOF(da, localC, &concs);
	checkPetscError(ierr,
		"PetscSolver0DHandler::computeDiagonalJacobian: "
		"DMDAVecRestoreKokkosOffsetViewDOF failed.");
}

} /* end namespace handler */
} /* end namespace solver */
} /* end namespace xolotl */
