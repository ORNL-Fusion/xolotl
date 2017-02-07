/**
 * An example of loading clusters
 */
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<Reactant.h>
#include<PSIClusterNetworkLoader.h>

using namespace std;

//! This operation prints proper usage instructions for this example 
void printUsage() {
	cout << "Usage:" << endl;
	cout << "\tclusterLoading <network file name>" << endl;
}

//! Main program
int main(int argc, char **argv) {

	// Local Declarations
	std::shared_ptr<ifstream> inputFileStream(new ifstream());
	xolotlCore::PSIClusterNetworkLoader clusterLoader = xolotlCore::PSIClusterNetworkLoader();

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


	return EXIT_SUCCESS;
}
