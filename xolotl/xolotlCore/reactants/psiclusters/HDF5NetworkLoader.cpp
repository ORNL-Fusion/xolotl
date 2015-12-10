#include "HDF5NetworkLoader.h"
#include <stdio.h>
#include <limits>
#include <algorithm>
#include <vector>
#include "HeCluster.h"
#include "VCluster.h"
#include "InterstitialCluster.h"
#include "HeVCluster.h"
#include "SuperCluster.h"
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

	// Apply sectional grouping
	applySectionalGrouping(network);

	return network;
}

void HDF5NetworkLoader::applySectionalGrouping(std::shared_ptr<PSIClusterReactionNetwork> network) {
	// Define the starting helium size
	int heMin = 31;
	// Define the width of the sections
	int sectionWidth = 2;

	// Get the HeV cluster map
	auto heVMap = network->getAll(heVType);

	// Create a temporary vector for the loop
	std::vector<PSICluster *> tempVector;

	// Initialize variables for the loop
	PSICluster * cluster;
	std::shared_ptr<SuperCluster> superCluster;
	static std::map<std::string, int> composition;
	int count = 0;
	double heSize = 0.0, vSize = 0.0, radius = 0.0, energy = 0.0;

	// Loop on the HeV clusters
	for (int i = 0; i < heVMap.size(); i++) {
		// Get the cluster and its composition
		cluster = (PSICluster *) heVMap.at(i);
		composition = cluster->getComposition();

		// Skip the clusters that are too small
		if (composition[heType] < heMin) continue;

		// Increment the counter
		count++;

		// Add this cluster to the temporary vector
		tempVector.push_back(cluster);
		heSize += (double) composition[heType];
		vSize += (double) composition[vType];
		radius += cluster->getReactionRadius();
		energy += cluster->getFormationEnergy();

		// Create the super cluster when the count is equal to
		// the width of the section
		if (count == sectionWidth) {
			// Average all values
			heSize = heSize / (double) sectionWidth;
			vSize = vSize / (double) sectionWidth;
			radius = radius / (double) sectionWidth;
			energy = energy / (double) sectionWidth;
			// Create the cluster
			superCluster = std::make_shared<SuperCluster>(heSize, vSize,
					sectionWidth, radius, energy, handlerRegistry);
			// Set the HeV vector
			superCluster->setHeVVector(tempVector);
			// Add this cluster to the network and clusters
			network->addSuper(superCluster);
			clusters.push_back(superCluster);

			// Reinitialize everything
			heSize = 0.0, vSize = 0.0, radius = 0.0, energy = 0.0;
			count = 0;
			tempVector.clear();
		}
	}

	// Gather the last clusters together if the count was not round
	if (count != 0) {
		// Average all values
		heSize = heSize / (double) count;
		vSize = vSize / (double) count;
		radius = radius / (double) count;
		energy = energy / (double) count;
		// Create the cluster
		superCluster = std::make_shared<SuperCluster>(heSize, vSize,
				count, radius, energy, handlerRegistry);
		// Set the HeV vector
		superCluster->setHeVVector(tempVector);
		// Add this cluster to the network and clusters
		network->addSuper(superCluster);
		clusters.push_back(superCluster);
	}

	// Initialize variables for the loop
	PSICluster * newCluster;
	// Loop on all the reactants to update the pairs vector with super clusters
	auto reactants = network->getAll();
	// Get the super cluster map
	auto superMap = network->getAll(superType);
	for (int i = 0; i < reactants->size(); i++) {
		// Get the cluster
		cluster = (PSICluster *) reactants->at(i);
		// Get their production and dissociation vectors
		auto react = cluster->reactingPairs;
		auto combi = cluster->combiningReactants;
		auto disso = cluster->dissociatingPairs;
		auto emi = cluster->emissionPairs;

		// Loop on its reacting pairs
		for (int l = 0; l < react.size(); l++) {
			// Test the first reactant
			if (react[l].first->getType() == heVType) {
				// Get its composition
				composition = react[l].first->getComposition();
				// Test its size
				if (composition[heType] >= heMin) {
					// It has to be replaced by a super cluster
					newCluster = findSuperCluster(superMap, composition);
					react[l].first = newCluster;
				}
			}

			// Test the second reactant
			if (react[l].second->getType() == heVType) {
				// Get its composition
				composition = react[l].second->getComposition();
				// Test its size
				if (composition[heType] >= heMin) {
					// It has to be replaced by a super cluster
					newCluster = findSuperCluster(superMap, composition);
					react[l].second = newCluster;
				}
			}
		}

		// Loop on its combining reactants
		for (int l = 0; l < combi.size(); l++) {
			// Test the combining reactant
			if (combi[l].combining->getType() == heVType) {
				// Get its composition
				composition = combi[l].combining->getComposition();
				// Test its size
				if (composition[heType] >= heMin) {
					// It has to be replaced by a super cluster
					newCluster = findSuperCluster(superMap, composition);
					combi[l].combining = newCluster;
				}
			}
		}

		// Loop on its dissociating pairs
		for (int l = 0; l < disso.size(); l++) {
			// Test the first reactant
			if (disso[l].first->getType() == heVType) {
				// Get its composition
				composition = disso[l].first->getComposition();
				// Test its size
				if (composition[heType] >= heMin) {
					// It has to be replaced by a super cluster
					newCluster = findSuperCluster(superMap, composition);
					disso[l].first = newCluster;
				}
			}

			// Test the second reactant
			if (disso[l].second->getType() == heVType) {
				// Get its composition
				composition = disso[l].second->getComposition();
				// Test its size
				if (composition[heType] >= heMin) {
					// It has to be replaced by a super cluster
					newCluster = findSuperCluster(superMap, composition);
					disso[l].second = newCluster;
				}
			}
		}

		// Loop on its emission pairs
		for (int l = 0; l < emi.size(); l++) {
			// Test the first reactant
			if (emi[l].first->getType() == heVType) {
				// Get its composition
				composition = emi[l].first->getComposition();
				// Test its size
				if (composition[heType] >= heMin) {
					// It has to be replaced by a super cluster
					newCluster = findSuperCluster(superMap, composition);
					emi[l].first = newCluster;
				}
			}

			// Test the second reactant
			if (emi[l].second->getType() == heVType) {
				// Get its composition
				composition = emi[l].second->getComposition();
				// Test its size
				if (composition[heType] >= heMin) {
					// It has to be replaced by a super cluster
					newCluster = findSuperCluster(superMap, composition);
					emi[l].second = newCluster;
				}
			}
		}

		// Set their production and dissociation vectors
		cluster->reactingPairs = react;
		cluster->combiningReactants = combi;
		cluster->dissociatingPairs = disso;
		cluster->emissionPairs = emi;
	}

	// Set the reaction network for each super reactant
	for (auto reactantsIt = superMap.begin(); reactantsIt != superMap.end();
			++reactantsIt) {
		(*reactantsIt)->setReactionNetwork(network);
	}

	// Remove HeV clusters bigger than heMin from the network
	// Loop on the HeV clusters
	for (int i = 0; i < heVMap.size(); i++) {
		// Get the cluster and its composition
		cluster = (PSICluster *) heVMap.at(i);
		composition = cluster->getComposition();

		// Skip the clusters that are too small
		if (composition[heType] < heMin) continue;

		// Reset the cluster
		cluster->reset();
	}

	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeNetwork();

	return;
}

PSICluster * HDF5NetworkLoader::findSuperCluster(std::vector<Reactant *> clusterList,
		std::map<std::string, int> comp) {
	// Initialize variables for the loop
	SuperCluster * superCluster;
	static std::map<std::string, int> tempComp;
	bool found = false;

	// Loop on the super cluster list
	for (int i = 0; i < clusterList.size(); i++) {
		// Get the super cluster
		superCluster = (SuperCluster *) clusterList.at(i);
		found = false;

		// Loop on its heVVector
		for (int l = 0; l < superCluster->heVVector.size(); l++) {
			// Get the composition and compare it
			tempComp = superCluster->heVVector[l]->getComposition();
			if (comp[heType] == tempComp[heType]
					&& comp[vType] == tempComp[vType]
					&& comp[iType] == tempComp[iType]) {
				// We got the right one!
				found = true;
				break;
			}
		}

		// If we found it we can return the super cluster
		if (found) break;
	}

	return (PSICluster *) superCluster;
}

void HDF5NetworkLoader::setFilename (const std::string& name) {
	fileName = name;

	return;
}

std::string HDF5NetworkLoader::getFilename () const {
	return fileName;
}

