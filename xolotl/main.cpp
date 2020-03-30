#include <interface.h>

//! Main program
int main(int argc, char **argv) {

	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Create an interface to control the solver
	XolotlInterface interface;

	// Initialize it
	interface.initializeXolotl(argc, argv);

	// Initialize the GB location
	interface.initGBLocation();

	// Run the solve
	interface.solveXolotl();

	// Finalize the run
	interface.finalizeXolotl();

	// Finalize MPI
	MPI_Finalize();

	return EXIT_SUCCESS;
}
