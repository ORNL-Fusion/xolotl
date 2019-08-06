#include <AlloyFaultedCluster.h>
#include <AlloyFrankCluster.h>
#include <AlloyIntCluster.h>
#include <AlloyPerfectCluster.h>
#include "AlloyClusterNetworkLoader.h"
#include <AlloySuperCluster.h>
#include <AlloyVacCluster.h>
#include <AlloyVoidCluster.h>
#include <xolotlPerf.h>
#include "xolotlCore/io/XFile.h"

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
		cluster = new AlloyVacCluster(numV, network, handlerRegistry);
	} else if (numI > 0) {
		cluster = new AlloyIntCluster(numI, network, handlerRegistry);
	} else if (numVoid > 0) {
		cluster = new AlloyVoidCluster(numVoid, network, handlerRegistry);
	} else if (numFaulted > 0) {
		cluster = new AlloyFaultedCluster(numFaulted, network, handlerRegistry);
	} else if (numFrank > 0) {
		cluster = new AlloyFrankCluster(numFrank, network, handlerRegistry);
	} else if (numPerfect > 0) {
		cluster = new AlloyPerfectCluster(numPerfect, network, handlerRegistry);
	}
	assert(cluster != nullptr);

	// TODO when we have widespread C++14 support, use std::make_unique
	// and construct unique ptr and object pointed to in one memory operation.
	return std::unique_ptr<AlloyCluster>(cluster);
}

std::unique_ptr<AlloyCluster> AlloyClusterNetworkLoader::createAlloySuperCluster(
		int nTot, int maxXe, ReactantType type,
		IReactionNetwork& network) const {
	// Create the cluster
	auto superCluster = new AlloySuperCluster(maxXe, nTot, type, network,
			handlerRegistry);

	// TODO when we have widespread C++14 support, use std::make_unique
	// and construct unique ptr and object pointed to in one memory operation.
	return std::unique_ptr<AlloyCluster>(superCluster);
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
	sizeMax = -1;
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
	sizeMax = -1;
	sectionWidth = 1;

	return;
}

std::unique_ptr<IReactionNetwork> AlloyClusterNetworkLoader::load(
		const IOptions& options) {
	// Get the dataset from the HDF5 files
	int normalSize = 0, superSize = 0;
	XFile networkFile(fileName);
	auto networkGroup = networkFile.getGroup<XFile::NetworkGroup>();
	assert(networkGroup);
	networkGroup->readNetworkSize(normalSize, superSize);

	// Initialization
	int numV = 0, numI = 0, numVoid = 0, numFaulted = 0, numFrank = 0,
			numPerfect = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// Prepare the network
	std::unique_ptr<AlloyClusterReactionNetwork> network(
			new AlloyClusterReactionNetwork(handlerRegistry));

	// Loop on the clusters
	for (int i = 0; i < normalSize + superSize; i++) {
		// Open the cluster group
		XFile::ClusterGroup clusterGroup(*networkGroup, i);

		if (i < normalSize) {
			// Normal cluster
			// Read the composition
			auto comp = clusterGroup.readCluster(formationEnergy,
					migrationEnergy, diffusionFactor);
			numV = comp[toCompIdx(Species::V)];
			numI = comp[toCompIdx(Species::I)];
			numFaulted = comp[toCompIdx(Species::Faulted)];
			numFrank = comp[toCompIdx(Species::Frank)];
			numVoid = comp[toCompIdx(Species::Void)];
			numPerfect = comp[toCompIdx(Species::Perfect)];

			// Create the cluster
			auto nextCluster = createAlloyCluster(numV, numI, numVoid,
					numFaulted, numFrank, numPerfect, *network);

			// Set the formation energy
			nextCluster->setFormationEnergy(formationEnergy);
			// Set the diffusion factor and migration energy
			nextCluster->setMigrationEnergy(migrationEnergy);
			nextCluster->setDiffusionFactor(diffusionFactor);

			// Save it in the network
			pushAlloyCluster(network, reactants, nextCluster);
		} else {
			// Super cluster
			int nTot = 0, maxAtom = 0;
			ReactantType type;
			clusterGroup.readAlloySuperCluster(nTot, maxAtom, type);

			// Create the cluster
			auto nextCluster = createAlloySuperCluster(nTot, maxAtom, type,
					*network);

			// Save it in the network
			pushAlloyCluster(network, reactants, nextCluster);
		}
	}

	// Ask reactants to update now that they are in network.
	for (IReactant& currReactant : reactants) {
		currReactant.updateFromNetwork();
	}

	// Set the reactions
	networkGroup->readReactions(*network);

	// Recompute Ids and network size
	network->reinitializeNetwork();

	// Need to use move() because return type uses smart pointer to base class,
	// not derived class that we created.
	// Some C++11 compilers accept it without the move, but apparently
	// that is not correct behavior until C++14.
	return std::move(network);
}

std::unique_ptr<IReactionNetwork> AlloyClusterNetworkLoader::generate(
		const IOptions &options) {
	// Initial declarations
	int maxV = options.getMaxV(), maxI = options.getMaxI();
	sizeMax = options.getMaxImpurity();
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
	for (int i = maxV + 1; i <= min(sizeMin - 1, sizeMax); ++i) {
		// Set the composition
		numVoid = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, 0, numVoid, 0, 0, 0, *network);

		// Save it in the network
		pushAlloyCluster(network, reactants, nextCluster);
	}

	// Generate the Faulted clusters
	for (int i = maxV + 1; i <= min(sizeMin - 1, sizeMax); ++i) {
		// Set the composition
		numFaulted = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, 0, 0, numFaulted, 0, 0,
				*network);

		// Save it in the network
		pushAlloyCluster(network, reactants, nextCluster);
	}

	// Generate the Frank clusters
	for (int i = maxI + 1; i <= min(sizeMin - 1, sizeMax); ++i) {
		// Set the composition
		numFrank = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, 0, 0, 0, numFrank, 0,
				*network);

		// Save it in the network
		pushAlloyCluster(network, reactants, nextCluster);
	}

	// Generate the Perfect clusters
	for (int i = maxI + 1; i <= min(sizeMin - 1, sizeMax); ++i) {
		// Set the composition
		numPerfect = i;
		// Create the cluster
		auto nextCluster = createAlloyCluster(0, 0, 0, 0, 0, numPerfect,
				*network);

		// Save it in the network
		pushAlloyCluster(network, reactants, nextCluster);
	}

	// Check if we want dummy reactions
	if (!dummyReactions) {
		// Apply sectional grouping
		applyGrouping(*network);
	}

	// Create the reactions
	network->createReactionConnectivity();

	return std::move(network);
}

void AlloyClusterNetworkLoader::applyGrouping(IReactionNetwork& network) const {

	// Decide here which types will undergo grouping
	std::vector<ReactantType> typeVec { ReactantType::Void,
			ReactantType::Faulted, ReactantType::Frank, ReactantType::Perfect };

	// Loop on them
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {
		auto currType = *tvIter;

		// Initialize variables for the loop
		std::unique_ptr<AlloySuperCluster> superCluster;
		int count = 0, superCount = 0, width = sectionWidth;
		int size = 0;

		// Loop on the xenon groups
		for (int k = sizeMin; k < sizeMax; k++) {

			// Increment the counter
			count++;

			// Track the size
			size = k;

			// Continue if we are not at the wanted width yet
			if (count < width && k < sizeMax - 1)
				continue;

			// Create the cluster
			auto rawSuperCluster = new AlloySuperCluster(size, count, currType,
					network, handlerRegistry);

//			std::cout << superCount << " " << count << " "
//					<< rawSuperCluster->getName() << std::endl;

			superCluster = std::unique_ptr<AlloySuperCluster>(rawSuperCluster);
			// Give the cluster to the network.
			network.add(std::move(superCluster));

			// Reinitialize everything
			size = 0;
			count = 0;
			superCount++;
//			width = std::max((int) (std::pow((double) superCount, 1.0) / 1.0),
//					sectionWidth);
//			width -= width % sectionWidth;
		}

		if (sizeMin < sizeMax) {
			// Group the last one alone
			auto rawSuperCluster = new AlloySuperCluster(sizeMax, 1, currType,
					network, handlerRegistry);

//		std::cout << superCount << " last " << rawSuperCluster->getName()
//				<< std::endl;

			superCluster = std::unique_ptr<AlloySuperCluster>(rawSuperCluster);
			// Give the cluster to the network.
			network.add(std::move(superCluster));
		}
	}

	// Recompute Ids and network size
	network.reinitializeNetwork();

	return;
}
