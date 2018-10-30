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

	// Print the retention
	interface.printRetention(solver);
	// Get the pointer to the vector
	auto retention = interface.getPointerRetention(solver);
	// Modify the vector
	if (retention->size() > 0)
		retention->at(3) = 1.0;
	// Print again
	interface.printRetention(solver);
	// Get the copy of the vector
	auto retentionCopy = interface.getCopyRetention(solver);
	// Modify the vector
	if (retentionCopy.size() > 0)
		retentionCopy[5] = 1.0;
	// Print again
	interface.printRetention(solver);

	// Finalize the run
	interface.finalizeXolotl(solver);

	// Finalize MPI
	MPI_Finalize();

	return EXIT_SUCCESS;
}
