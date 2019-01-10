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
	int xs, xm, Mx, ys, ym, My, zs, zm, Mz;
	interface.getLocalCoordinates(solver, xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	std::cout << "The local coordinates are: " << xs << " " << xm << " "
			<< Mx << "; " << ys << " " << ym << " "
			<< My << "; " << zs << " " << zm << " "
			<< Mz << std::endl;
	auto localRate = interface.getLocalXeRate(solver);

	for (auto k = 0; k < max(zm, 1); k++) {
		for (int j = 0; j < max(ym, 1); j++) {
			for (int i = 0; i < max(xm, 1); i++) {
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
