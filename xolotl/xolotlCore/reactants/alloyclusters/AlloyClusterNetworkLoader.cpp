#include "AlloyClusterNetworkLoader.h"
#include <VacCluster.h>
#include <IntCluster.h>
#include <VoidCluster.h>
#include <FaultedCluster.h>
#include <FrankCluster.h>
#include <PerfectCluster.h>
#include <AlloySuperCluster.h>
#include <xolotlPerf.h>

using namespace xolotlCore;

std::unique_ptr<AlloyCluster> AlloyClusterNetworkLoader::createAlloyCluster(
		int numV, int numI, int numVoid, int numFaulted, int numFrank,
		int numPerfect, IReactionNetwork& network) {
	// Local Declarations
	AlloyCluster* cluster;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if (numV > 0) {
		cluster = new VacCluster(numV, network, handlerRegistry);
	} else if (numI > 0) {
		cluster = new IntCluster(numI, network, handlerRegistry);
	} else if (numVoid > 0) {
		cluster = new VoidCluster(numVoid, network, handlerRegistry);
	} else if (numFaulted > 0) {
		cluster = new FaultedCluster(numFaulted, network, handlerRegistry);
	} else if (numFrank > 0) {
		cluster = new FrankCluster(numFrank, network, handlerRegistry);
	} else if (numPerfect > 0) {
		cluster = new PerfectCluster(numPerfect, network, handlerRegistry);
	}
	assert(cluster != nullptr);

	// TODO when we have widespread C++14 support, use std::make_unique
	// and construct unique ptr and object pointed to in one memory operation.
	return std::unique_ptr<AlloyCluster>(cluster);
}

void AlloyClusterNetworkLoader::pushAlloyCluster(
		std::unique_ptr<AlloyClusterReactionNetwork> & network,
		std::vector<std::reference_wrapper<Reactant> > & reactants,
		std::unique_ptr<AlloyCluster> & cluster) {
	// Check if we want dummy reactions
	if (dummyReactions) {
		// Create a dummy cluster (Reactant) from the existing cluster
		auto dummyCluster = std::unique_ptr<Reactant>(new Reactant(*cluster));
		// Save access to it so we can trigger updates after
		// we add all to the network.
		reactants.emplace_back(*dummyCluster);

		// Give the cluster to the network
		network->add(std::move(dummyCluster));
	} else {
		// Save access to it so we can trigger updates after
		// we add all to the network.
		reactants.emplace_back(*cluster);

		// Give the cluster to the network
		network->add(std::move(cluster));
	}

	return;
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

std::unique_ptr<IReactionNetwork> AlloyClusterNetworkLoader::load(
		const IOptions& options) {
	// Prepare the network
	// Once we have C++14 support, use std::make_unique.
	std::unique_ptr<AlloyClusterReactionNetwork> network(
			new AlloyClusterReactionNetwork(handlerRegistry));

	return network;
}

std::unique_ptr<IReactionNetwork> AlloyClusterNetworkLoader::generate(
		const IOptions &options) {
	// Initial declarations
	int maxV = options.getMaxV(), maxI = options.getMaxI(), maxBig =
			options.getMaxImpurity();
	int numV = 0, numI = 0, numVoid = 0, numFaulted = 0, numFrank = 0,
			numPerfect = 0;
	// Once we have C++14 support, use std::make_unique.
	std::unique_ptr<AlloyClusterReactionNetwork> network(
			new AlloyClusterReactionNetwork(handlerRegistry));
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// Generate the Vacancy clusters
	for (int i = 1; i <= maxV; ++i) {
		// Set the composition
		numV = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(numV, 0, 0, 0, 0, 0, *network);

		// Save it in the network
		pushAlloyCluster(network, reactants, nextCluster);
	}

	// Generate the Interstitial clusters
	for (int i = 1; i <= maxI; ++i) {
		// Set the composition
		numI = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, numI, 0, 0, 0, 0, *network);

		// Save it in the network
		pushAlloyCluster(network, reactants, nextCluster);
	}

	// Generate the Void clusters
	for (int i = maxV + 1; i <= maxBig; ++i) {
		// Set the composition
		numVoid = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, 0, numVoid, 0, 0, 0, *network);

		// Save it in the network
		pushAlloyCluster(network, reactants, nextCluster);
	}

	// Generate the Faulted clusters
	for (int i = maxV + 1; i <= maxBig; ++i) {
		// Set the composition
		numFaulted = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, 0, 0, numFaulted, 0, 0,
				*network);

		// Save it in the network
		pushAlloyCluster(network, reactants, nextCluster);
	}

	// Generate the Frank clusters
	for (int i = maxI + 1; i <= maxBig; ++i) {
		// Set the composition
		numFrank = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, 0, 0, 0, numFrank, 0,
				*network);

		// Save it in the network
		pushAlloyCluster(network, reactants, nextCluster);
	}

	// Generate the Perfect clusters
	for (int i = maxI + 1; i <= maxBig; ++i) {
		// Set the composition
		numPerfect = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, 0, 0, 0, 0, numPerfect,
				*network);

		// Save it in the network
		pushAlloyCluster(network, reactants, nextCluster);
	}

	// Ask reactants to update now that they are in network.
	for (IReactant& currReactant : reactants) {
		currReactant.updateFromNetwork();
	}

	// Create the reactions
	network->createReactionConnectivity();

	// Check if we want dummy reactions
	if (!dummyReactions) {
		// Apply sectional grouping
		applyGrouping(*network);
	}

	return network;
}

void AlloyClusterNetworkLoader::applyGrouping(IReactionNetwork& network) {

	// Decide here which types will undergo grouping
	std::vector<ReactantType> typeVec { ReactantType::Void,
			ReactantType::Faulted, ReactantType::Frank, ReactantType::Perfect };

	// Loop on them
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {
		auto currType = *tvIter;

		// Get the biggest size
		int sizeMax = network.getMaxClusterSize(currType);

		// Create a temporary vector for the loop
		std::vector<AlloyCluster *> tempVector;

		// Initialize variables for the loop
		AlloyCluster * cluster;
		std::unique_ptr<AlloySuperCluster> superCluster;
		int count = 0, superCount = 0, width = sectionWidth;
		double size = 0.0, radius = 0.0, energy = 0.0;

		// Map to know which cluster is in which group
		std::map<int, int> clusterGroupMap;
		// Map to know which super cluster gathers which group
		std::map<int, AlloySuperCluster *> superGroupMap;

		// Loop on the xenon groups
		for (int k = sizeMin; k <= sizeMax; k++) {
			// Get the corresponding cluster
			cluster = (AlloyCluster *) network.get(toSpecies(currType), k);

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
			if (count < width && k < sizeMax)
				continue;

			// Average all values
			size = size / (double) count;
			radius = radius / (double) count;
			energy = energy / (double) count;

			// Create the cluster
			auto rawSuperCluster = new AlloySuperCluster(size, count, count,
					radius, energy, currType, network, handlerRegistry);
			auto superCluster = std::unique_ptr<AlloySuperCluster>(
					rawSuperCluster);
			// Save access to the cluster so we can trigger updates
			// after we give it to the network.
			auto& scref = *superCluster;
			// Set the HeV vector
			scref.setAtomVector(tempVector);
			// Give the cluster to the network.
			network.add(std::move(superCluster));

			// Keep the information of the group
			superGroupMap[superCount] = rawSuperCluster;

//			std::cout << superCount << " " << count << " "
//					<< rawSuperCluster->getName() << std::endl;

			// Reinitialize everything
			size = 0.0, radius = 0.0, energy = 0.0;
			count = 0;
			tempVector.clear();
			superCount++;
//			width = std::max((int) (std::pow((double) superCount, 1.0) / 1.0),
//					sectionWidth);
//			width -= width % sectionWidth;
		}

		// Initialize variables for the loop
		AlloySuperCluster * newCluster;
		int clusterSize = 0;

		// Tell each reactant to update the pairs vector with super clusters
		for (IReactant& currReactant : network.getAll()) {

			auto& cluster = static_cast<AlloyCluster&>(currReactant);
			// Get their production and dissociation vectors
			auto react = cluster.reactingPairs;
			auto combi = cluster.combiningReactants;
			auto disso = cluster.dissociatingPairs;
			auto emi = cluster.emissionPairs;

			// Loop on its reacting pairs
			for (int l = 0; l < react.size(); l++) {
				// Test the first reactant
				if (react[l].first->getType() == currType) {
					// Get its size
					clusterSize = react[l].first->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster =
								superGroupMap[clusterGroupMap[clusterSize]];
						react[l].first = newCluster;
						react[l].firstDistance = newCluster->getDistance(
								clusterSize);
					}
				}

				// Test the second reactant
				if (react[l].second->getType() == currType) {
					// Get its size
					clusterSize = react[l].second->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster =
								superGroupMap[clusterGroupMap[clusterSize]];
						react[l].second = newCluster;
						react[l].secondDistance = newCluster->getDistance(
								clusterSize);
					}
				}
			}

			// Loop on its combining reactants
			for (int l = 0; l < combi.size(); l++) {
				// Test the combining reactant
				if (combi[l].combining->getType() == currType) {
					// Get its size
					clusterSize = combi[l].combining->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster =
								superGroupMap[clusterGroupMap[clusterSize]];
						combi[l].combining = newCluster;
						combi[l].distance = newCluster->getDistance(
								clusterSize);
					}
				}
			}

			// Loop on its dissociating pairs
			for (int l = 0; l < disso.size(); l++) {
				// Test the first reactant
				if (disso[l].first->getType() == currType) {
					// Get its size
					clusterSize = disso[l].first->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster =
								superGroupMap[clusterGroupMap[clusterSize]];
						disso[l].first = newCluster;
						disso[l].firstDistance = newCluster->getDistance(
								clusterSize);
					}
				}

				// Test the second reactant
				if (disso[l].second->getType() == currType) {
					// Get its size
					clusterSize = disso[l].second->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster =
								superGroupMap[clusterGroupMap[clusterSize]];
						disso[l].second = newCluster;
						disso[l].secondDistance = newCluster->getDistance(
								clusterSize);
					}
				}
			}

			// Loop on its emission pairs
			for (int l = 0; l < emi.size(); l++) {
				// Test the first reactant
				if (emi[l].first->getType() == currType) {
					// Get its size
					clusterSize = emi[l].first->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster =
								superGroupMap[clusterGroupMap[clusterSize]];
						emi[l].first = newCluster;
						emi[l].firstDistance = newCluster->getDistance(
								clusterSize);
					}
				}

				// Test the second reactant
				if (emi[l].second->getType() == currType) {
					// Get its size
					clusterSize = emi[l].second->getSize();
					// Test its size
					if (clusterSize >= sizeMin) {
						// It has to be replaced by a super cluster
						newCluster =
								superGroupMap[clusterGroupMap[clusterSize]];
						emi[l].second = newCluster;
						emi[l].secondDistance = newCluster->getDistance(
								clusterSize);
					}
				}
			}

			// Set their production and dissociation vectors
			cluster.reactingPairs = react;
			cluster.combiningReactants = combi;
			cluster.dissociatingPairs = disso;
			cluster.emissionPairs = emi;
		}
	}

	// Set the reaction network for each type of super reactant
	std::vector<ReactantType> superTypeVec { ReactantType::VoidSuper,
			ReactantType::FaultedSuper, ReactantType::FrankSuper,
			ReactantType::PerfectSuper };

	// Loop on them
	for (auto tvIter = superTypeVec.begin(); tvIter != superTypeVec.end();
			++tvIter) {
		auto currType = *tvIter;
		for (auto const& superMapItem : network.getAll(currType)) {
			auto& currCluster =
					static_cast<AlloySuperCluster&>(*(superMapItem.second));
			currCluster.updateFromNetwork();
		}
	}

	// Remove the clusters bigger than sizeMin from the network
	// Loop on the simple types
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {
		auto currType = *tvIter;
		std::vector<std::reference_wrapper<IReactant> > doomedReactants;
		for (auto const& currMapItem : network.getAll(currType)) {
			auto& currCluster = currMapItem.second;

			// Get the cluster's size.
			int clusterSize = currCluster->getSize();

			// Check if the cluster is too large.
			if (clusterSize >= sizeMin) {
				// The cluster is too large.  Add it to the ones we will remove.
				doomedReactants.push_back(*currCluster);
			}
		}
		network.removeReactants(doomedReactants);
	}

	// Recompute Ids and network size
	network.reinitializeNetwork();

	return;
}
