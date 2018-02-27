#include "HDF5NetworkLoader.h"
#include <stdio.h>
#include <limits>
#include <algorithm>
#include <vector>
#include "PSIClusterReactionNetwork.h"
#include <xolotlPerf.h>
#include <HDF5Utils.h>

using namespace xolotlCore;

std::unique_ptr<IReactionNetwork> HDF5NetworkLoader::load(
		const IOptions& options) {
	// Get the dataset from the HDF5 files
	auto networkVector = xolotlCore::HDF5Utils::readNetwork(fileName);

	// Initialization
	int numHe = 0, numV = 0, numI = 0, numW = 0, numD = 0, numT = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector < std::reference_wrapper<Reactant> > reactants;

	// Prepare the network
	std::unique_ptr < PSIClusterReactionNetwork
			> network(new PSIClusterReactionNetwork(handlerRegistry));

	// Loop on the networkVector
	for (auto lineIt = networkVector.begin(); lineIt != networkVector.end();
			lineIt++) {
		// Composition of the cluster
		numHe = (int) (*lineIt)[0];
		numD = (int) (*lineIt)[1];
		numT = (int) (*lineIt)[2];
		numV = (int) (*lineIt)[3];
		numI = (int) (*lineIt)[4];

		// If the cluster is big enough to be grouped
		if (numV >= vMin && numHe > 0) {
			// Created the coordinates and store them
			auto pair = std::make_pair(numHe, numV);
			heVList.emplace(pair);
		} else {
			// Create the cluster
			auto nextCluster = createPSICluster(numHe, numD, numT, numV, numI,
					*network);

			// Energies
			formationEnergy = (*lineIt)[5];
			migrationEnergy = (*lineIt)[6];
			diffusionFactor = (*lineIt)[7];

			// Set the formation energy
			nextCluster->setFormationEnergy(formationEnergy);
			// Set the diffusion factor and migration energy
			nextCluster->setMigrationEnergy(migrationEnergy);
			nextCluster->setDiffusionFactor(diffusionFactor);

			// Save it in the network
			pushPSICluster(network, reactants, nextCluster);
		}
	}

	// Ask reactants to update now that they are in network.
	for (IReactant& currReactant : reactants) {
		currReactant.updateFromNetwork();
	}

	// Check if we want dummy reactions
	if (!dummyReactions) {
		// Apply sectional grouping
		applySectionalGrouping (*network);
	}

	// Create the reactions
	network->createReactionConnectivity();

	// Recompute Ids and network size
	network->reinitializeNetwork();

	// Need to use move() because return type uses smart pointer to base class,
	// not derived class that we created.
	// Some C++11 compilers accept it without the move, but apparently
	// that is not correct behavior until C++14.
	return std::move(network);
}
