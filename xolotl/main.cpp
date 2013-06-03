/**
 * Main.c, currently only able to load clusters
 */
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<Reactant.h>
#include<PSIClusterNetworkLoader.h>
#include<PetscSolver.h>

using namespace std;

//! This operation prints the start message
void printStartMessage() {
	cout << "Starting Xolotl Plasma-Surface Interactions Simulator" << endl;
	// TODO! Print copyright message
	// TODO! Print date and time
}

//! This operation prints proper usage instructions for Xolotl
void printUsage() {
	cout << "Usage:" << endl;
	cout << "\txolotl <network file name>" << endl;
}

//! Main program
int main(int argc, char **argv) {

	// Local Declarations
	std::shared_ptr<ifstream> inputFileStream(new ifstream());
	xolotlCore::PSIClusterNetworkLoader clusterLoader =
			xolotlCore::PSIClusterNetworkLoader();
	xolotlSolver::PetscSolver solver;

	// Print the start message
	printStartMessage();

	// Check the arguments
	if (argc < 2) {
		cout << "Insufficient input provided! Aborting!" << std::endl;
		printUsage();
		return EXIT_FAILURE;
	}

	// Load the file from the input argument
	inputFileStream->open(argv[1]);

	try {
		// Load the clusters
		clusterLoader.setInputstream(inputFileStream);
		shared_ptr<xolotlCore::ReactionNetwork> network = clusterLoader.load();
		// Print some otherwise useless information to show it worked.
		std::map<std::string, std::string> props = *(network->properties);
		cout << "Loaded input file, " << argv[1] << endl;
		cout << "Found:" << endl;
		cout << "\t" << props["numHeClusters"]
				<< " helium clusters with max cluster size "
				<< props["maxHeClusterSize"] << endl;
		cout << "\t" << props["numVClusters"]
				<< " vacancy clusters with max cluster size "
				<< props["maxVClusterSize"] << endl;
		cout << "\t" << props["numIClusters"]
				<< " interstitial clusters with max cluster size "
				<< props["maxIClusterSize"] << endl;
		cout << "\t" << props["numMixedClusters"]
				<< " mixed-species clusters with max cluster size " << endl;

		// Setup and run the solver
		solver.setCommandLineOptions(argc,argv);
		solver.initialize();
		solver.solve();
		solver.finalize();

	} catch (std::string error) {
		std::cout << error << std::endl;
		std::cout << "Aborting." << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
