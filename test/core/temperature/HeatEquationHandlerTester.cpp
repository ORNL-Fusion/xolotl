#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <mpi.h>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/Constants.h>
#include <xolotl/core/temperature/HeatEquationHandler.h>
#include <xolotl/options/Options.h>
#include <xolotl/test/Util.h>
#include <xolotl/test/config.h>

using namespace std;
using namespace xolotl;
using namespace core;
using namespace temperature;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * This suite is responsible for testing the HeatEquationHandler.
 */
BOOST_AUTO_TEST_SUITE(HeatEquationHandler_testSuite)

/**
 * Method checking the initialization of the off-diagonal and diagonal part of
 * the Jacobian, and the compute temperature methods.
 */
BOOST_AUTO_TEST_CASE(checkHeat1D)
{
	MPI_Init(NULL, NULL);
	// Set the DOF
	const int dof = 9;

	// Create the heat handler
	auto heatHandler = HeatEquationHandler(5.0e-12, 1000.0, 1);
	heatHandler.setHeatCoefficient(tungstenHeatCoefficient);
	heatHandler.setHeatConductivity(tungstenHeatConductivity);

	// Check the initial temperatures
	BOOST_REQUIRE_CLOSE(
		heatHandler.getTemperature({0.0, 0.0, 0.0}, 0.0), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		heatHandler.getTemperature({1.0, 0.0, 0.0}, 0.0), 1000.0, 0.01);

	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double)l);
	}

	// Set a time
	double time = 0.5;

	// Initialize it
	heatHandler.initialize(dof);
	heatHandler.updateSurfacePosition(0, grid);

	// The size parameter in the x direction
	double hx = 1.0;

	// The arrays of concentration
	test::DOFView concentration("concentration", 3, dof + 1);
	test::DOFView newConcentration("newConcentration", 3, dof + 1);

	// Initialize their values
	Kokkos::parallel_for(
		Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {3, dof + 1}),
		KOKKOS_LAMBDA(int i, int n) {
			auto id = static_cast<double>(i * (dof + 1) + n);
			concentration(i, n) = id * id;
		});

	// Get the offset for the grid point in the middle
	// Supposing the 3 grid points are laid-out as follow:
	// 0 | 1 | 2
	auto concOffset = subview(concentration, 1, Kokkos::ALL);
	auto updatedConcOffset = subview(newConcentration, 1, Kokkos::ALL);

	// Fill the concVector with the pointer to the middle, left, and right grid
	// points
	using ConcSubView = Kokkos::View<const double*>;
	Kokkos::Array<ConcSubView, 3> concVector;
	concVector[0] = concOffset; // middle
	concVector[1] = subview(concentration, 0, Kokkos::ALL); // left
	concVector[2] = subview(concentration, 2, Kokkos::ALL); // right

	// Compute the heat equation at this grid point
	heatHandler.computeTemperature(
		time, concVector.data(), updatedConcOffset, hx, hx, hx);

	// Check the new values of updatedConcOffset
	auto updatedConcOffsetMirror =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[9], 7.5004e+15, 0.01);

	// Set the temperature in the handler
	heatHandler.setTemperature(concOffset);
	// Check the updated temperature
	plsm::SpaceVector<double, 3> pos{1.0, 0.0, 0.0};
	BOOST_REQUIRE_CLOSE(heatHandler.getTemperature(pos, 1.0), 361.0, 0.01);

	// Initialize the indices and values to set in the Jacobian
	IdType indices[1];
	double val[3];
	// Get the pointer on them for the compute diffusion method
	IdType* indicesPointer = &indices[0];
	double* valPointer = &val[0];

	Kokkos::Array<ConcSubView::host_mirror_type, 3> hConcVec;
	const double* hConcPtrVec[3];
	int id = 0;
	for (auto&& xId : {1, 0, 2}) {
		concVector[id] = subview(concentration, xId, Kokkos::ALL);
		hConcVec[id] = create_mirror_view(concVector[id]);
		deep_copy(hConcVec[id], concVector[id]);
		hConcPtrVec[id] = hConcVec[id].data();
		++id;
	}

	// Compute the partial derivatives for the heat equation a the grid point
	heatHandler.computePartialsForTemperature(
		time, hConcPtrVec, valPointer, indicesPointer, hx, hx, hx);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 9);

	// Check the values
	BOOST_REQUIRE_CLOSE(val[0], -135508639961461, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 80171695638351, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 50744437570026, 0.01);
}

BOOST_AUTO_TEST_CASE(checkHeat2D)
{
	// Set the DOF
	const int dof = 9;

	// Create the heat handler
	auto heatHandler = HeatEquationHandler(5.0e-12, 1000.0, 2);
	heatHandler.setHeatCoefficient(tungstenHeatCoefficient);
	heatHandler.setHeatConductivity(tungstenHeatConductivity);

	// Check the initial temperatures
	BOOST_REQUIRE_CLOSE(
		heatHandler.getTemperature({0.0, 0.0, 0.0}, 0.0), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		heatHandler.getTemperature({1.0, 0.0, 0.0}, 0.0), 1000.0, 0.01);

	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double)l);
	}

	// Set a time
	double time = 0.5;

	// Initialize it
	heatHandler.initialize(dof);
	heatHandler.updateSurfacePosition(0, grid);

	// The step size in the x direction
	double hx = 1.0;
	// The size parameter in the y direction
	double sy = 1.0;

	// The arrays of concentration
	test::DOFView concentration("concentration", 9, dof + 1);
	test::DOFView newConcentration("newConcentration", 9, dof + 1);

	// Initialize their values
	Kokkos::parallel_for(
		Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {9, dof + 1}),
		KOKKOS_LAMBDA(int i, int n) {
			auto id = static_cast<double>(i * (dof + 1) + n);
			concentration(i, n) = id * id;
		});

	// Get the offset for the grid point in the middle
	// Supposing the 9 grid points are laid-out as follow:
	// 6 | 7 | 8
	// 3 | 4 | 5
	// 0 | 1 | 2
	auto concOffset = subview(concentration, 4, Kokkos::ALL);
	auto updatedConcOffset = subview(newConcentration, 4, Kokkos::ALL);

	// Fill the concVector with the pointer to the middle, left, right, bottom,
	// and top grid points
	using ConcSubView = Kokkos::View<const double*>;
	Kokkos::Array<ConcSubView, 5> concVector;
	concVector[0] = concOffset; // middle
	concVector[1] = subview(concentration, 3, Kokkos::ALL); // left
	concVector[2] = subview(concentration, 5, Kokkos::ALL); // right
	concVector[3] = subview(concentration, 1, Kokkos::ALL); // bottom
	concVector[4] = subview(concentration, 7, Kokkos::ALL); // top

	// Compute the heat equation at this grid point
	heatHandler.computeTemperature(
		time, concVector.data(), updatedConcOffset, hx, hx, hx, sy, 1);

	// Check the new values of updatedConcOffset
	auto updatedConcOffsetMirror =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[9], 1.367e+17, 0.01);

	// Set the temperature in the handler
	heatHandler.setTemperature(concOffset);
	// Check the updated temperature
	plsm::SpaceVector<double, 3> pos{1.0, 0.0, 0.0};
	BOOST_REQUIRE_CLOSE(heatHandler.getTemperature(pos, 1.0), 2401, 0.01);

	// Initialize the indices and values to set in the Jacobian
	IdType indices[1];
	double val[5];
	// Get the pointer on them for the compute diffusion method
	IdType* indicesPointer = &indices[0];
	double* valPointer = &val[0];

	Kokkos::Array<ConcSubView::host_mirror_type, 5> hConcVec;
	const double* hConcPtrVec[5];
	int id = 0;
	for (auto&& xId : {4, 3, 5, 1, 7}) {
		concVector[id] = subview(concentration, xId, Kokkos::ALL);
		hConcVec[id] = create_mirror_view(concVector[id]);
		deep_copy(hConcVec[id], concVector[id]);
		hConcPtrVec[id] = hConcVec[id].data();
		++id;
	}

	// Compute the partial derivatives for the heat equation a the grid point
	heatHandler.computePartialsForTemperature(
		time, hConcPtrVec, valPointer, indicesPointer, hx, hx, hx, sy, 1);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 9);

	// Check the values
	BOOST_REQUIRE_CLOSE(val[0], -109442547402517, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 34626084704135, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 8451742066439, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 21538913385287, 0.01);
	BOOST_REQUIRE_CLOSE(val[4], 21538913385287, 0.01);
}

BOOST_AUTO_TEST_CASE(checkHeat3D)
{
	// Set the DOF
	const int dof = 9;

	// Create the heat handler
	auto heatHandler = HeatEquationHandler(5.0e-12, 1000.0, 3);
	heatHandler.setHeatCoefficient(tungstenHeatCoefficient);
	heatHandler.setHeatConductivity(tungstenHeatConductivity);

	// Check the initial temperatures
	BOOST_REQUIRE_CLOSE(
		heatHandler.getTemperature({0.0, 0.0, 0.0}, 0.0), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		heatHandler.getTemperature({1.0, 0.0, 0.0}, 0.0), 1000.0, 0.01);

	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double)l);
	}

	// Set a time
	double time = 0.5;

	// Initialize it
	heatHandler.initialize(dof);
	heatHandler.updateSurfacePosition(0, grid);

	// The step size in the x direction
	double hx = 1.0;
	// The size parameter in the y direction
	double sy = 1.0;
	// The size parameter in the z direction
	double sz = 1.0;

	// The arrays of concentration
	test::DOFView concentration("concentration", 27, dof + 1);
	test::DOFView newConcentration("newConcentration", 27, dof + 1);

	// Initialize their values
	Kokkos::parallel_for(
		Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {27, dof + 1}),
		KOKKOS_LAMBDA(int i, int n) {
			auto id = static_cast<double>(i * (dof + 1) + n);
			concentration(i, n) = id * id / 10.0;
		});

	// Get the offset for the grid point in the middle
	// Supposing the 27 grid points are laid-out as follow (a cube!):
	// 6 | 7 | 8    15 | 16 | 17    24 | 25 | 26
	// 3 | 4 | 5    12 | 13 | 14    21 | 22 | 23
	// 0 | 1 | 2    9  | 10 | 11    18 | 19 | 20
	//   front         middle           back
	auto concOffset = subview(concentration, 13, Kokkos::ALL);
	auto updatedConcOffset = subview(newConcentration, 13, Kokkos::ALL);

	// Fill the concVector with the pointer to the middle, left, right, bottom,
	// top, front, and back grid points
	using ConcSubView = Kokkos::View<const double*>;
	Kokkos::Array<ConcSubView, 7> concVector;
	concVector[0] = concOffset; // middle
	concVector[1] = subview(concentration, 12, Kokkos::ALL); // left
	concVector[2] = subview(concentration, 14, Kokkos::ALL); // right
	concVector[3] = subview(concentration, 10, Kokkos::ALL); // bottom
	concVector[4] = subview(concentration, 16, Kokkos::ALL); // top
	concVector[5] = subview(concentration, 4, Kokkos::ALL); // front
	concVector[6] = subview(concentration, 22, Kokkos::ALL); // back

	// Compute the heat equation at this grid point
	heatHandler.computeTemperature(
		time, concVector.data(), updatedConcOffset, hx, hx, hx, sy, 1, sz, 1);

	// Check the new values of updatedConcOffset
	auto updatedConcOffsetMirror =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[9], 1.24397e+17, 0.01);

	// Set the temperature in the handler
	heatHandler.setTemperature(concOffset);
	// Check the updated temperature
	plsm::SpaceVector<double, 3> pos{1.0, 0.0, 0.0};
	BOOST_REQUIRE_CLOSE(heatHandler.getTemperature(pos, 1.0), 1932.1, 0.01);

	// Initialize the indices and values to set in the Jacobian
	IdType indices[1];
	double val[7];
	// Get the pointer on them for the compute diffusion method
	IdType* indicesPointer = &indices[0];
	double* valPointer = &val[0];

	Kokkos::Array<ConcSubView::host_mirror_type, 7> hConcVec;
	const double* hConcPtrVec[7];
	int id = 0;
	for (auto&& xId : {13, 12, 14, 10, 16, 4, 22}) {
		concVector[id] = subview(concentration, xId, Kokkos::ALL);
		hConcVec[id] = create_mirror_view(concVector[id]);
		deep_copy(hConcVec[id], concVector[id]);
		hConcPtrVec[id] = hConcVec[id].data();
		++id;
	}

	// Compute the partial derivatives for the heat equation a the grid point
	heatHandler.computePartialsForTemperature(time, hConcPtrVec, valPointer,
		indicesPointer, hx, hx, hx, sy, 1, sz, 1);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 9);

	// Check the values
	BOOST_REQUIRE_CLOSE(val[0], -164490411906709, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 29368449647960, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 25268952000254, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 27318700824107, 0.01);
	BOOST_REQUIRE_CLOSE(val[4], 27318700824107, 0.01);
	BOOST_REQUIRE_CLOSE(val[5], 27318700824107, 0.01);
	BOOST_REQUIRE_CLOSE(val[6], 27318700824107, 0.01);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
