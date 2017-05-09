#include "AlloyClusterNetworkLoader.h"
#include <AlloyClusterReactionNetwork.h>
#include <VacCluster.h>
#include <IntCluster.h>
#include <VoidCluster.h>
#include <FaultedCluster.h>
#include <FrankCluster.h>
#include <PerfectCluster.h>
#include <AlloySuperCluster.h>
#include <HDF5Utils.h>
#include <xolotlPerf.h>

using namespace xolotlCore;

std::shared_ptr<AlloyCluster> AlloyClusterNetworkLoader::createAlloyCluster(
		int numV, int numI, int numVoid, int numFaulted, int numFrank,
		int numPerfect) {
	// Local Declarations
	std::shared_ptr<AlloyCluster> cluster;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if (numV > 0) {
		cluster = std::make_shared<VacCluster>(numV, handlerRegistry);
	} else if (numI > 0) {
		cluster = std::make_shared<IntCluster>(numI, handlerRegistry);
	} else if (numVoid > 0) {
		cluster = std::make_shared<VoidCluster>(numVoid, handlerRegistry);
	} else if (numFaulted > 0) {
		cluster = std::make_shared<FaultedCluster>(numFaulted, handlerRegistry);
	} else if (numFrank > 0) {
		cluster = std::make_shared<FrankCluster>(numFrank, handlerRegistry);
	} else if (numPerfect > 0) {
		cluster = std::make_shared<PerfectCluster>(numPerfect, handlerRegistry);
	}

	return cluster;
}

AlloyClusterNetworkLoader::AlloyClusterNetworkLoader(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = nullptr;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	sizeMin = 1000000;
	sectionWidth = 1;

	return;
}

AlloyClusterNetworkLoader::AlloyClusterNetworkLoader(
		const std::shared_ptr<std::istream> stream,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = stream;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	sizeMin = 1000000;
	sectionWidth = 1;

	return;
}

std::shared_ptr<IReactionNetwork> AlloyClusterNetworkLoader::load() {
	// Prepare the network
	std::shared_ptr<AlloyClusterReactionNetwork> network = std::make_shared<
			AlloyClusterReactionNetwork>(handlerRegistry);

	return network;
}

std::shared_ptr<IReactionNetwork> AlloyClusterNetworkLoader::generate(
		IOptions &options) {
	// Initial declarations
	int maxV = options.getMaxV(), maxI = options.getMaxI(), maxBig =
			options.getMaxImpurity();
	int numV = 0, numI = 0, numVoid = 0, numFaulted = 0, numFrank = 0,
			numPerfect = 0;
	std::shared_ptr<AlloyClusterReactionNetwork> network = std::make_shared<
			AlloyClusterReactionNetwork>(handlerRegistry);
	std::vector<std::shared_ptr<Reactant> > reactants;

	// THE DIFFUSION PARAMETERS HAVE TO BE DEFINED HERE

	// Generate the Vacancy clusters
	for (int i = 1; i <= maxV; ++i) {
		// Set the composition
		numV = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(numV, 0, 0, 0, 0, 0);

		// THE DIFFUSION PARAMETERS HAVE TO BE SET HERE

		// Add the cluster to the network
		network->add(nextCluster);
		// Add it to the list so that we can set the network later
		reactants.push_back(nextCluster);
	}

	// Generate the Interstitial clusters
	for (int i = 1; i <= maxI; ++i) {
		// Set the composition
		numI = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, numI, 0, 0, 0, 0);

		// THE DIFFUSION PARAMETERS HAVE TO BE SET HERE

		// Add the cluster to the network
		network->add(nextCluster);
		// Add it to the list so that we can set the network later
		reactants.push_back(nextCluster);
	}

	// Generate the Void clusters
	for (int i = maxV + 1; i <= maxBig; ++i) {
		// Set the composition
		numVoid = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, 0, numVoid, 0, 0, 0);

		// THE DIFFUSION PARAMETERS HAVE TO BE SET HERE

		// Add the cluster to the network
		network->add(nextCluster);
		// Add it to the list so that we can set the network later
		reactants.push_back(nextCluster);
	}

	// Generate the Faulted clusters
	for (int i = maxV + 1; i <= maxBig; ++i) {
		// Set the composition
		numFaulted = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, 0, 0, numFaulted, 0, 0);

		// THE DIFFUSION PARAMETERS HAVE TO BE SET HERE

		// Add the cluster to the network
		network->add(nextCluster);
		// Add it to the list so that we can set the network later
		reactants.push_back(nextCluster);
	}

	// Generate the Frank clusters
	for (int i = maxI + 1; i <= maxBig; ++i) {
		// Set the composition
		numFrank = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, 0, 0, 0, numFrank, 0);

		// THE DIFFUSION PARAMETERS HAVE TO BE SET HERE

		// Add the cluster to the network
		network->add(nextCluster);
		// Add it to the list so that we can set the network later
		reactants.push_back(nextCluster);
	}

	// Generate the Perfect clusters
	for (int i = maxI + 1; i <= maxBig; ++i) {
		// Set the composition
		numPerfect = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, 0, 0, 0, 0, numPerfect);

		// THE DIFFUSION PARAMETERS HAVE TO BE SET HERE

		// Add the cluster to the network
		network->add(nextCluster);
		// Add it to the list so that we can set the network later
		reactants.push_back(nextCluster);
	}

	// Set the network for all of the reactants. This MUST be done manually.
	for (auto currCluster : reactants) {
		currCluster->setReactionNetwork(network);
	}

	// Create the reactions
	network->createReactionConnectivity();

	// Check if we want dummy reactions
	if (!dummyReactions) {
		// Apply sectional grouping
		applyGrouping(network);
	}

	return network;
}

void AlloyClusterNetworkLoader::applyGrouping(
		std::shared_ptr<IReactionNetwork> network) {
	//Initial declarations
	std::vector<IReactant*> doomedReactants;

	// Decide here which types will undergo grouping
	std::vector<string> typeVec = { voidType, faultedType, frankType,
			perfectType };

	// Loop on them
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {
		string typeName = *tvIter;

		// Get all the reactants of this type
		auto allTypeReactants = network->getAll(typeName);

		// Create a temporary vector for the loop
		std::vector<AlloyCluster *> tempVector;

		// Initialize variables for the loop
		AlloyCluster * cluster;
		std::shared_ptr<AlloySuperCluster> superCluster;
		int count = 0, superCount = 0, width = sectionWidth;
		double size = 0.0, radius = 0.0, energy = 0.0;

		// Map to know which cluster is in which group
		std::map<int, int> clusterGroupMap;
		// Map to know which super cluster gathers which group
		std::map<int, AlloySuperCluster *> superGroupMap;

		// Loop on the xenon groups
		for (int k = sizeMin; k <= allTypeReactants.size(); k++) {
			// Get the corresponding cluster
			cluster = (AlloyCluster *) network->get(typeName, k);

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
			if (count < width && k < allTypeReactants.size())
				continue;

			// Average all values
			size = size / (double) count;
			radius = radius / (double) count;
			energy = energy / (double) count;
			// Create the cluster
			superCluster = std::make_shared<AlloySuperCluster>(size, count,
					count, radius, energy, typeName, handlerRegistry);
			// Set the HeV vector
			superCluster->setAtomVector(tempVector);
			// Add this cluster to the network and clusters
			network->addSuper(superCluster);
			// Keep the information of the group
			superGroupMap[superCount] = superCluster.get();

//			std::cout << superCount << " " << count << " "
//					<< superCluster->getName() << std::endl;

			// Reinitialize everything
			size = 0.0, radius = 0.0, energy = 0.0;
			count = 0;
			tempVector.clear();
			superCount++;
			width = std::max((int) std::pow((double) superCount, 3.0),
					sectionWidth);
			width -= width % sectionWidth;
		}

		// Initialize variables for the loop
		AlloySuperCluster * newCluster;
		int clusterSize = 0;
		// Loop on all the reactants to update the pairs vector with super clusters
		auto reactants = network->getAll();
		for (int i = 0; i < reactants->size(); i++) {
			// Get the cluster
			cluster = (AlloyCluster *) reactants->at(i);
			// Get their production and dissociation vectors
			auto react = cluster->reactingPairs;
			auto combi = cluster->combiningReactants;
			auto disso = cluster->dissociatingPairs;
			auto emi = cluster->emissionPairs;

			// Loop on its reacting pairs
			for (int l = 0; l < react.size(); l++) {
				// Test the first reactant
				if (react[l].first->getType() == typeName) {
					// Get its size
					clusterSize = react[l].first->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster = superGroupMap[clusterGroupMap[clusterSize]];
						react[l].first = newCluster;
						react[l].firstDistance = newCluster->getDistance(clusterSize);
					}
				}

				// Test the second reactant
				if (react[l].second->getType() == typeName) {
					// Get its size
					clusterSize = react[l].second->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster = superGroupMap[clusterGroupMap[clusterSize]];
						react[l].second = newCluster;
						react[l].secondDistance = newCluster->getDistance(clusterSize);
					}
				}
			}

			// Loop on its combining reactants
			for (int l = 0; l < combi.size(); l++) {
				// Test the combining reactant
				if (combi[l].combining->getType() == typeName) {
					// Get its size
					clusterSize = combi[l].combining->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster = superGroupMap[clusterGroupMap[clusterSize]];
						combi[l].combining = newCluster;
						combi[l].distance = newCluster->getDistance(clusterSize);
					}
				}
			}

			// Loop on its dissociating pairs
			for (int l = 0; l < disso.size(); l++) {
				// Test the first reactant
				if (disso[l].first->getType() == typeName) {
					// Get its size
					clusterSize = disso[l].first->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster = superGroupMap[clusterGroupMap[clusterSize]];
						disso[l].first = newCluster;
						disso[l].firstDistance = newCluster->getDistance(clusterSize);
					}
				}

				// Test the second reactant
				if (disso[l].second->getType() == typeName) {
					// Get its size
					clusterSize = disso[l].second->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster = superGroupMap[clusterGroupMap[clusterSize]];
						disso[l].second = newCluster;
						disso[l].secondDistance = newCluster->getDistance(clusterSize);
					}
				}
			}

			// Loop on its emission pairs
			for (int l = 0; l < emi.size(); l++) {
				// Test the first reactant
				if (emi[l].first->getType() == typeName) {
					// Get its size
					clusterSize = emi[l].first->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster = superGroupMap[clusterGroupMap[clusterSize]];
						emi[l].first = newCluster;
						emi[l].firstDistance = newCluster->getDistance(clusterSize);
					}
				}

				// Test the second reactant
				if (emi[l].second->getType() == typeName) {
					// Get its size
					clusterSize = emi[l].second->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster = superGroupMap[clusterGroupMap[clusterSize]];
						emi[l].second = newCluster;
						emi[l].secondDistance = newCluster->getDistance(clusterSize);
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
		auto superMap = network->getAll(AlloySuperType);
		// Set the reaction network for each super reactant
		for (auto currCluster : superMap) {
			currCluster->setReactionNetwork(network);
		}

		// Remove clusters bigger than sizeMin from the network
		// Loop on the clusters
		for (auto currCluster : allTypeReactants) {

			// Get the cluster's size.
			clusterSize = currCluster->getSize();

			// Check if the cluster is too large.
			if (clusterSize >= sizeMin) {
				// The cluster is too large.  Add it to the ones we will remove.
				doomedReactants.push_back(currCluster);
			}
		}
	}
	network->removeReactants(doomedReactants);

	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeNetwork();

	return;
}
