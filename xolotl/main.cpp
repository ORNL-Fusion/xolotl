#include <xolotl/interface/Interface.h>

//! Main program
int main(int argc, char **argv) {

	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Initialize kokkos
	Kokkos::initialize();

	// Create an interface to control the solver
	xolotl::interface::XolotlInterface interface = xolotl::interface::XolotlInterface();

	// Initialize it
	interface.initializeXolotl(argc, argv);

	// Initialize the GB location
	interface.initGBLocation();

	// Run the solve
	interface.solveXolotl();

	// Finalize the run
	interface.finalizeXolotl();

	// Finalize kokkos
	Kokkos::finalize();

	// Finalize MPI
	MPI_Finalize();

	return EXIT_SUCCESS;
}
