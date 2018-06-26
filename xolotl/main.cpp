#include <interface.h>

//! Main program
int main(int argc, char **argv) {
	// Create an interface to control the solver
	XolotlInterface interface;

	// Initialize it
	auto solver = interface.initializeXolotl(argc, argv);
	// Run the solve
	interface.solveXolotl(solver);
	// Finalize the run
	interface.finalizeXolotl(solver);

	return EXIT_SUCCESS;
}
