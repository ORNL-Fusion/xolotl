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
#include <HandlerRegistryFactory.h>
#include <HDF5Utils.h>

using namespace xolotlCore;

std::shared_ptr<PSIClusterReactionNetwork> HDF5NetworkLoader::load() {
	// Get the dataset from the HDF5 files
	auto networkVector = xolotlCore::HDF5Utils::readNetwork(fileName);

	// Initialization
	int numHe = 0, numV = 0, numI = 0;
	double heBindingE = 0.0, vBindingE = 0.0, iBindingE = 0.0, migrationEnergy =
			0.0;
	double diffusionFactor = 0.0;
	std::vector<double> bindingEnergies;
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

		// Binding energies
		heBindingE = (*lineIt)[3];
		vBindingE = (*lineIt)[4];
		iBindingE = (*lineIt)[5];
		migrationEnergy = (*lineIt)[6];
		diffusionFactor = (*lineIt)[7];

		// Create the binding energies array and set it
		bindingEnergies.clear();
		bindingEnergies.push_back(heBindingE);
		bindingEnergies.push_back(vBindingE);
		bindingEnergies.push_back(iBindingE);
		nextCluster->setBindingEnergies(bindingEnergies);
		// Set the diffusion factor and migration energy
		nextCluster->setMigrationEnergy(migrationEnergy);
		nextCluster->setDiffusionFactor(diffusionFactor);
		// Add the cluster to the network
		network->add(nextCluster);
		// Add it to the list so that we can set the network later
		reactants.push_back(nextCluster);
	}

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

