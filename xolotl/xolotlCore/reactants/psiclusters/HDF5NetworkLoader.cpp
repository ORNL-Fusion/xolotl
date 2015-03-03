#include "HDF5NetworkLoader.h"
#include <stdio.h>
#include <limits>
#include <algorithm>
#include <vector>
#include "HeCluster.h"
#include "VCluster.h"
#include "InterstitialCluster.h"
#include "HeVCluster.h"
// #include "HeInterstitialCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <xolotlPerf.h>
#include <HDF5Utils.h>

using namespace xolotlCore;

std::shared_ptr<PSIClusterReactionNetwork> HDF5NetworkLoader::load() {
	// Get the dataset from the HDF5 files
	auto networkVector = xolotlCore::HDF5Utils::readNetwork(fileName);

	// Initialization
	int numHe = 0, numV = 0, numI = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::shared_ptr<Reactant> > reactants;

	// Prepare the network
	std::shared_ptr<PSIClusterReactionNetwork> network = std::make_shared
			< PSIClusterReactionNetwork > (handlerRegistry);

	// Loop on the networkVector
	for (auto lineIt = networkVector.begin(); lineIt != networkVector.end();
			lineIt++) {
		// Composition of the cluster
		numHe = (int) (*lineIt)[0];
		numV = (int) (*lineIt)[1];
		numI = (int) (*lineIt)[2];
		// Create the cluster
		auto nextCluster = createCluster(numHe, numV, numI);

		// Energies
		formationEnergy = (*lineIt)[3];
		migrationEnergy = (*lineIt)[4];
		diffusionFactor = (*lineIt)[5];

		// Set the formation energy
		nextCluster->setFormationEnergy(formationEnergy);
		// Set the diffusion factor and migration energy
		nextCluster->setMigrationEnergy(migrationEnergy);
		nextCluster->setDiffusionFactor(diffusionFactor);
		// Add the cluster to the network
		network->add(nextCluster);
		// Add it to the list so that we can set the network later
		reactants.push_back(nextCluster);
	}

	// Set the reaction network for each reactant
	for (auto reactantsIt = reactants.begin(); reactantsIt != reactants.end();
			++reactantsIt) {
		(*reactantsIt)->setReactionNetwork(network);
	}

	return network;
}

void HDF5NetworkLoader::setFilename (std::string name) {
	fileName = name;

	return;
}

std::string HDF5NetworkLoader::getFilename () const {
	return fileName;
}

