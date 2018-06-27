#include <interface.h>

//! Main program
int main(int argc, char **argv) {

	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Create an interface to control the solver
	XolotlInterface interface;

	// Initialize it
	auto solver = interface.initializeXolotl(argc, argv);
	// Run the solve
	interface.solveXolotl(solver);
	// Finalize the run
	interface.finalizeXolotl(solver);

	// Finalize MPI
	MPI_Finalize();

	return EXIT_SUCCESS;
}
