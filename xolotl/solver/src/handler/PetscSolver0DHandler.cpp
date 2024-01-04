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
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (size > 1) {
		throw std::runtime_error("\nYou are trying to run a 0D simulation in "
								 "parallel, this is not possible!");
	}

	PetscCallVoid(
		DMDACreate1d(xolotlComm, DM_BOUNDARY_NONE, 1, dof + 1, 0, NULL, &da));
	PetscCallVoid(DMSetFromOptions(da));
	PetscCallVoid(DMSetUp(da));
}

void
PetscSolver0DHandler::initializeSolverContext(DM& da, TS& ts)
{
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

	// Tell the network the number of grid points on this process
	network.setGridSize(1);

	// Get the diagonal fill
	auto nPartials = network.getDiagonalFill(dfill);

	// Preallocate matrix
	Mat J;
	PetscCallVoid(TSGetRHSJacobian(ts, &J, nullptr, nullptr, nullptr));
	auto [rows, cols] = convertToCoordinateListPair(dof, dfill);
	// handling temperature (FIXME)
	rows.push_back(dof);
	cols.push_back(dof);
	++nPartials;
	//
	PetscCallVoid(
		MatSetPreallocationCOO(J, rows.size(), rows.data(), cols.data()));

	// Initialize the arrays for the reaction partial derivatives
	vals = Kokkos::View<double*>("solverPartials", nPartials);

	// Set the size of the partial derivatives vectors
	reactingPartialsForCluster.resize(dof, 0.0);

	// Initialize the flux handler
	fluxHandler->initializeFluxHandler(network, 0, grid);
}

void
PetscSolver0DHandler::initializeConcentration(
	DM& da, Vec& C, DM& oldDA, Vec& oldC)
{
	// Initialize the last temperature
	temperature.push_back(0.0);

	// Pointer for the concentration vector
	PetscScalar** concentrations = nullptr;
	PetscCallVoid(DMDAVecGetArrayDOF(da, C, &concentrations));

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
	if (not hasConcentrations) {
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
	PetscCallVoid(DMDAVecRestoreArrayDOF(da, C, &concentrations));

	return;
}

std::vector<std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
PetscSolver0DHandler::getConcVector(DM& da, Vec& C)
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
	PetscCallContinue(DMDAVecRestoreArrayDOFRead(da, C, &concentrations));

	return toReturn;
}

void
PetscSolver0DHandler::setConcVector(DM& da, Vec& C,
	std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>&
		concVector)
{
	// Pointer for the concentration vector
	PetscScalar* gridPointSolution = nullptr;
	PetscScalar** concentrations = nullptr;
	PetscCallVoid(DMDAVecGetArrayDOF(da, C, &concentrations));

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
	PetscCallVoid(DMDAVecRestoreArrayDOF(da, C, &concentrations));

	return;
}

void
PetscSolver0DHandler::updateConcentration(
	TS& ts, Vec& localC, Vec& F, PetscReal ftime)
{
	// Get the local data vector from PETSc
	DM da;
	PetscCallVoid(TSGetDM(ts, &da));

	// Pointers to the PETSc arrays that start at the beginning of the
	// local array
	PetscOffsetView<const PetscScalar**> concs;
	PetscCallVoid(DMDAVecGetKokkosOffsetViewDOF(da, localC, &concs));
	PetscOffsetView<PetscScalar**> updatedConcs;
	PetscCallVoid(DMDAVecGetKokkosOffsetViewDOFWrite(da, F, &updatedConcs));

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
	PetscCallVoid(DMDAVecRestoreKokkosOffsetViewDOF(da, localC, &concs));
	PetscCallVoid(DMDAVecRestoreKokkosOffsetViewDOFWrite(da, F, &updatedConcs));
}

void
PetscSolver0DHandler::computeJacobian(
	TS& ts, Vec& localC, Mat& J, PetscReal ftime)
{
	// Get the distributed array
	DM da;
	PetscCallVoid(TSGetDM(ts, &da));

	// Get pointers to vector data
	PetscOffsetView<const PetscScalar**> concs;
	PetscCallVoid(DMDAVecGetKokkosOffsetViewDOF(da, localC, &concs));

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
