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

	// Check the local Xe rate at a few grid points
	auto localRate = interface.getLocalXeRate(solver);

	for (auto k = 0; k < localRate->at(0)[0].size(); k++) {
		for (int j = 0; j < localRate->at(0).size(); j++) {
			for (int i = 0; i < localRate->size(); i++) {
				std::cout << localRate->at(i)[j][k] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	// Finalize the run
	interface.finalizeXolotl(solver);

	// Finalize MPI
	MPI_Finalize();

	return EXIT_SUCCESS;
}
