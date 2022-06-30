#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <mpi.h>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/Constants.h>
#include <xolotl/core/temperature/HeatEquationHandler.h>
#include <xolotl/options/Options.h>
#include <xolotl/test/config.h>

using namespace std;
using namespace xolotl;
using namespace core;
using namespace temperature;

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

	// Create ofill
	network::IReactionNetwork::SparseFillMap ofill;
	// Create dfill
	network::IReactionNetwork::SparseFillMap dfill;

	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double)l);
	}

	// Set a time
	double time = 0.5;

	// Initialize it
	heatHandler.initializeTemperature(dof, ofill, dfill);

	// Check that the temperature "diffusion" is well set
	BOOST_REQUIRE_EQUAL(ofill[9][0], 9);
	BOOST_REQUIRE_EQUAL(dfill[9][0], 9);

	// The size parameter in the x direction
	double hx = 1.0;

	// The arrays of concentration
	double concentration[3 * (dof + 1)];
	double newConcentration[3 * (dof + 1)];

	// Initialize their values
	for (int i = 0; i < 3 * (dof + 1); i++) {
		concentration[i] = (double)i * i;
		newConcentration[i] = 0.0;
	}

	// Get pointers
	double* conc = &concentration[0];
	double* updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	// Supposing the 3 grid points are laid-out as follow:
	// 0 | 1 | 2
	double* concOffset = conc + (dof + 1);
	double* updatedConcOffset = updatedConc + (dof + 1);

	// Fill the concVector with the pointer to the middle, left, and right grid
	// points
	double* concVector[3]{};
	concVector[0] = concOffset; // middle
	concVector[1] = conc; // left
	concVector[2] = conc + 2 * (dof + 1); // right

	// Compute the heat equation at this grid point
	heatHandler.computeTemperature(
		time, concVector, updatedConcOffset, hx, hx, hx);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[9], 1.367e+16, 0.01);

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

	// Compute the partial derivatives for the heat equation a the grid point
	heatHandler.computePartialsForTemperature(
		time, concVector, valPointer, indicesPointer, hx, hx, hx);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 9);

	// Check the values
	BOOST_REQUIRE_CLOSE(val[0], -1.367e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 6.835e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 6.835e+13, 0.01);
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

	// Create ofill
	network::IReactionNetwork::SparseFillMap ofill;
	// Create dfill
	network::IReactionNetwork::SparseFillMap dfill;

	// Initialize it
	heatHandler.initializeTemperature(dof, ofill, dfill);

	// Check that the temperature "diffusion" is well set
	BOOST_REQUIRE_EQUAL(ofill[9][0], 9);
	BOOST_REQUIRE_EQUAL(dfill[9][0], 9);

	// The step size in the x direction
	double hx = 1.0;
	// The size parameter in the y direction
	double sy = 1.0;

	// The arrays of concentration
	double concentration[9 * (dof + 1)];
	double newConcentration[9 * (dof + 1)];

	// Initialize their values
	for (int i = 0; i < 9 * (dof + 1); i++) {
		concentration[i] = (double)i * i;
		newConcentration[i] = 0.0;
	}

	// Get pointers
	double* conc = &concentration[0];
	double* updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	// Supposing the 9 grid points are laid-out as follow:
	// 6 | 7 | 8
	// 3 | 4 | 5
	// 0 | 1 | 2
	double* concOffset = conc + 4 * (dof + 1);
	double* updatedConcOffset = updatedConc + 4 * (dof + 1);

	// Fill the concVector with the pointer to the middle, left, right, bottom,
	// and top grid points
	double* concVector[5]{};
	concVector[0] = concOffset; // middle
	concVector[1] = conc + 3 * (dof + 1); // left
	concVector[2] = conc + 5 * (dof + 1); // right
	concVector[3] = conc + 1 * (dof + 1); // bottom
	concVector[4] = conc + 7 * (dof + 1); // top

	// Compute the heat equation at this grid point
	heatHandler.computeTemperature(
		time, concVector, updatedConcOffset, hx, hx, hx, sy, 1);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[9], 1.367e+17, 0.01);

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

	// Compute the partial derivatives for the heat equation a the grid point
	heatHandler.computePartialsForTemperature(
		time, concVector, valPointer, indicesPointer, hx, hx, hx, sy, 1);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 9);

	// Check the values
	BOOST_REQUIRE_CLOSE(val[0], -2.734e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 6.835e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 6.835e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 6.835e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[4], 6.835e+13, 0.01);
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

	// Create ofill
	network::IReactionNetwork::SparseFillMap ofill;
	// Create dfill
	network::IReactionNetwork::SparseFillMap dfill;

	// Initialize it
	heatHandler.initializeTemperature(dof, ofill, dfill);

	// Check that the temperature "diffusion" is well set
	BOOST_REQUIRE_EQUAL(ofill[9][0], 9);
	BOOST_REQUIRE_EQUAL(dfill[9][0], 9);

	// The step size in the x direction
	double hx = 1.0;
	// The size parameter in the y direction
	double sy = 1.0;
	// The size parameter in the z direction
	double sz = 1.0;

	// The arrays of concentration
	double concentration[27 * (dof + 1)];
	double newConcentration[27 * (dof + 1)];

	// Initialize their values
	for (int i = 0; i < 27 * (dof + 1); i++) {
		concentration[i] = (double)i * i / 10.0;
		newConcentration[i] = 0.0;
	}

	// Get pointers
	double* conc = &concentration[0];
	double* updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	// Supposing the 27 grid points are laid-out as follow (a cube!):
	// 6 | 7 | 8    15 | 16 | 17    24 | 25 | 26
	// 3 | 4 | 5    12 | 13 | 14    21 | 22 | 23
	// 0 | 1 | 2    9  | 10 | 11    18 | 19 | 20
	//   front         middle           back
	double* concOffset = conc + 13 * (dof + 1);
	double* updatedConcOffset = updatedConc + 13 * (dof + 1);

	// Fill the concVector with the pointer to the middle, left, right, bottom,
	// top, front, and back grid points
	double* concVector[7]{};
	concVector[0] = concOffset; // middle
	concVector[1] = conc + 12 * (dof + 1); // left
	concVector[2] = conc + 14 * (dof + 1); // right
	concVector[3] = conc + 10 * (dof + 1); // bottom
	concVector[4] = conc + 16 * (dof + 1); // top
	concVector[5] = conc + 4 * (dof + 1); // front
	concVector[6] = conc + 22 * (dof + 1); // back

	// Compute the heat equation at this grid point
	heatHandler.computeTemperature(
		time, concVector, updatedConcOffset, hx, hx, hx, sy, 1, sz, 1);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[9], 1.24397e+17, 0.01);

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

	// Compute the partial derivatives for the heat equation a the grid point
	heatHandler.computePartialsForTemperature(
		time, concVector, valPointer, indicesPointer, hx, hx, hx, sy, 1, sz, 1);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 9);

	// Check the values
	BOOST_REQUIRE_CLOSE(val[0], -4.101e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 6.835e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 6.835e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 6.835e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[4], 6.835e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[5], 6.835e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[6], 6.835e+13, 0.01);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
