#include "NEClusterNetworkLoader.h"
#include <NEClusterReactionNetwork.h>
#include <XeCluster.h>
#include <NESuperCluster.h>
#include <HDF5Utils.h>
#include <xolotlPerf.h>

using namespace xolotlCore;

///**
// * This operation converts a string to a double, taking in to account the fact
// * that the input file may contain keys such as "infinite."
// *
// * @param inString the string to be converted
// * @return the string as a double
// */
//static inline double convertStrToDouble(const std::string& inString) {
//	return (inString.compare("infinite") == 0) ?
//			std::numeric_limits<double>::infinity() :
//			strtod(inString.c_str(), NULL);
//}

std::shared_ptr<NECluster> NEClusterNetworkLoader::createNECluster(int numXe,
		int numV, int numI) {
	// Local Declarations
	std::shared_ptr<NECluster> cluster;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if (numXe > 0) {
		// Create a new XeVCluster
		cluster = std::make_shared<XeCluster>(numXe, handlerRegistry);
	}

	return cluster;
}

NEClusterNetworkLoader::NEClusterNetworkLoader(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = nullptr;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	xeMin = 1000000;
	sectionWidth = 1;

	return;
}

NEClusterNetworkLoader::NEClusterNetworkLoader(
		const std::shared_ptr<std::istream> stream,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = stream;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	xeMin = 1000000;
	sectionWidth = 1;

	return;
}

std::shared_ptr<IReactionNetwork> NEClusterNetworkLoader::load() {
	// Get the dataset from the HDF5 files
	auto networkVector = xolotlCore::HDF5Utils::readNetwork(fileName);

	// Initialization
	int numXe = 0, numV = 0, numI = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::shared_ptr<Reactant> > reactants;

	// Prepare the network
	std::shared_ptr<NEClusterReactionNetwork> network = std::make_shared<
			NEClusterReactionNetwork>(handlerRegistry);

	// Loop on the networkVector
	for (auto lineIt = networkVector.begin(); lineIt != networkVector.end();
			lineIt++) {

		// Composition of the cluster
		numXe = (int) (*lineIt)[0];
		numV = (int) (*lineIt)[1];
		numI = (int) (*lineIt)[2];
		// Create the cluster
		auto nextCluster = createNECluster(numXe, numV, numI);

		// Energies
		formationEnergy = (*lineIt)[3];
		migrationEnergy = (*lineIt)[4];
		diffusionFactor = (*lineIt)[5];

		// Set the formation energy
		nextCluster->setFormationEnergy(formationEnergy);
		// Set the diffusion factor and migration energy
		nextCluster->setMigrationEnergy(migrationEnergy);
		nextCluster->setDiffusionFactor(diffusionFactor);

		// Check if we want dummy reactions
		if (dummyReactions) {
			// Create a dummy cluster (Reactant) from the existing cluster
			auto dummyCluster = std::static_pointer_cast<Reactant>(
					nextCluster->Reactant::clone());
			// Add the cluster to the network
			network->add(dummyCluster);
			// Add it to the list so that we can set the network later
			reactants.push_back(dummyCluster);
		} else {
			// Add the cluster to the network
			network->add(nextCluster);
			// Add it to the list so that we can set the network later
			reactants.push_back(nextCluster);
		}
	}

	// Set the reaction network for each reactant
	for (auto currCluster : reactants) {
		currCluster->setReactionNetwork(network);
	}

	// Check if we want dummy reactions
	if (!dummyReactions) {
		// Apply grouping
		applyGrouping(network);
	}

	return network;
}

void NEClusterNetworkLoader::applyGrouping(
		std::shared_ptr<IReactionNetwork> network) {
	// Get the xenon cluster map
	auto xeMap = network->getAll(xeType);

	// Create a temporary vector for the loop
	std::vector<NECluster *> tempVector;

	// Initialize variables for the loop
	NECluster * cluster;
	std::shared_ptr<NESuperCluster> superCluster;
	int count = 0, superCount = 0, width = sectionWidth;
	double size = 0.0, radius = 0.0, energy = 0.0;

	// Map to know which cluster is in which group
	std::map<int, int> clusterGroupMap;
	// Map to know which super cluster gathers which group
	std::map<int, NESuperCluster *> superGroupMap;

	// Loop on the xenon groups
	for (int k = xeMin; k <= xeMap.size(); k++) {
		// Get the corresponding cluster
		cluster = (NECluster *) network->get(xeType, k);

		// Verify if the cluster exists
		if (!cluster)
			continue;

		// Increment the counter
		count++;

		// Add this cluster to the temporary vector
		tempVector.push_back(cluster);
		size += (double) k;
		radius += cluster->getReactionRadius();
		energy += cluster->getFormationEnergy();

		// Save in which group it is
		clusterGroupMap[k] = superCount;

		// Check if there were clusters in this group
		if (count < width && k < xeMap.size())
			continue;

		// Average all values
		size = size / (double) count;
		radius = radius / (double) count;
		energy = energy / (double) count;
		// Create the cluster
		superCluster = std::make_shared<NESuperCluster>(size, count, count,
				radius, energy, handlerRegistry);
		// Set the HeV vector
		superCluster->setXeVector(tempVector);
		// Add this cluster to the network and clusters
		network->addSuper(superCluster);
		// Keep the information of the group
		superGroupMap[superCount] = superCluster.get();

//		std::cout << superCount << " " << count << " "
//				<< superCluster->getName() << std::endl;

		// Reinitialize everything
		size = 0.0, radius = 0.0, energy = 0.0;
		count = 0;
		tempVector.clear();
		superCount++;
		width = max(sectionWidth * (int) (superCount / 10), sectionWidth);
	}

	// Initialize variables for the loop
	NESuperCluster * newCluster;
	int nXe = 0;
	// Loop on all the reactants to update the pairs vector with super clusters
	auto reactants = network->getAll();
	for (int i = 0; i < reactants->size(); i++) {
		// Get the cluster
		cluster = (NECluster *) reactants->at(i);
		// Get their production and dissociation vectors
		auto react = cluster->reactingPairs;
		auto combi = cluster->combiningReactants;
		auto disso = cluster->dissociatingPairs;
		auto emi = cluster->emissionPairs;

		// Loop on its reacting pairs
		for (int l = 0; l < react.size(); l++) {
			// Test the first reactant
			if (react[l].first->getType() == xeType) {
				// Get its size
				nXe = react[l].first->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					newCluster = superGroupMap[clusterGroupMap[nXe]];
					react[l].first = newCluster;
					react[l].firstDistance = newCluster->getDistance(size);
				}
			}

			// Test the second reactant
			if (react[l].second->getType() == xeType) {
				// Get its size
				nXe = react[l].second->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					newCluster = superGroupMap[clusterGroupMap[nXe]];
					react[l].second = newCluster;
					react[l].secondDistance = newCluster->getDistance(nXe);
				}
			}
		}

		// Loop on its combining reactants
		for (int l = 0; l < combi.size(); l++) {
			// Test the combining reactant
			if (combi[l].combining->getType() == xeType) {
				// Get its size
				nXe = combi[l].combining->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					newCluster = superGroupMap[clusterGroupMap[nXe]];
					combi[l].combining = newCluster;
					combi[l].distance = newCluster->getDistance(nXe);
				}
			}
		}

		// Loop on its dissociating pairs
		for (int l = 0; l < disso.size(); l++) {
			// Test the first reactant
			if (disso[l].first->getType() == xeType) {
				// Get its size
				nXe = disso[l].first->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					newCluster = superGroupMap[clusterGroupMap[nXe]];
					disso[l].first = newCluster;
					disso[l].firstDistance = newCluster->getDistance(nXe);
				}
			}

			// Test the second reactant
			if (disso[l].second->getType() == xeType) {
				// Get its size
				nXe = disso[l].second->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					newCluster = superGroupMap[clusterGroupMap[nXe]];
					disso[l].second = newCluster;
					disso[l].secondDistance = newCluster->getDistance(nXe);
				}
			}
		}

		// Loop on its emission pairs
		for (int l = 0; l < emi.size(); l++) {
			// Test the first reactant
			if (emi[l].first->getType() == xeType) {
				// Get its size
				nXe = emi[l].first->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					newCluster = superGroupMap[clusterGroupMap[nXe]];
					emi[l].first = newCluster;
					emi[l].firstDistance = newCluster->getDistance(nXe);
				}
			}

			// Test the second reactant
			if (emi[l].second->getType() == xeType) {
				// Get its size
				nXe = emi[l].second->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					newCluster = superGroupMap[clusterGroupMap[nXe]];
					emi[l].second = newCluster;
					emi[l].secondDistance = newCluster->getDistance(nXe);
				}
			}
		}

		// Set their production and dissociation vectors
		cluster->reactingPairs = react;
		cluster->combiningReactants = combi;
		cluster->dissociatingPairs = disso;
		cluster->emissionPairs = emi;
	}

	// Get the super cluster map
	auto superMap = network->getAll(NESuperType);
	// Set the reaction network for each super reactant
	for (auto currCluster : superMap) {
		currCluster->setReactionNetwork(network);
	}

	// Remove Xe clusters bigger than xeMin from the network
	// Loop on the Xe clusters
	std::vector<IReactant*> doomedReactants;
	for (auto currCluster : xeMap) {

		// Get the cluster's size.
		nXe = currCluster->getSize();

		// Check if the cluster is too large.
		if (nXe >= xeMin) {
			// The cluster is too large.  Add it to the ones we will remove.
			doomedReactants.push_back(currCluster);
		}
	}
	network->removeReactants(doomedReactants);

	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeNetwork();

	return;
}
