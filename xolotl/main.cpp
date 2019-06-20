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

//	// Get the concentration vector at the end
//	auto concVector = interface.getConcVector();
//
//	// Loop on it
//	std::cout << concVector.size() << std::endl;
//	std::cout << concVector[0].size() << std::endl;
//	std::cout << concVector[0][0].size() << std::endl;
//	for (int k = 0; k < concVector.size(); k++) {
//		for (int j = 0; j < concVector[0].size(); j++) {
//			for (int i = 0; i < concVector[0][0].size(); i++) {
//				std::cout << i << " " << j << " " << k << " "
//						<< concVector[k][j][i].size() << " "
//						<< concVector[k][j][i][0].first << " "
//						<< concVector[k][j][i][0].second << std::endl;
//			}
//		}
//	}

	// Finalize the run
	interface.finalizeXolotl();

	// Finalize MPI
	MPI_Finalize();

	return EXIT_SUCCESS;
}
