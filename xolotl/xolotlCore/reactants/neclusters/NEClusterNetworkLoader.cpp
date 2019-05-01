#include <fstream>
#include <functional>
#include <cassert>
#include "NEClusterNetworkLoader.h"
#include <NEClusterReactionNetwork.h>
#include <NEXeCluster.h>
#include <NESuperCluster.h>
#include <MPIUtils.h>
#include <xolotlPerf.h>
#include "xolotlCore/io/XFile.h"

namespace xolotlCore {

std::unique_ptr<NECluster> NEClusterNetworkLoader::createNECluster(int numXe,
		int numV, int numI, IReactionNetwork& network) const {

	// Local Declarations
	NECluster* cluster = nullptr;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if (numXe > 0) {
		// Create a new XeVCluster
		cluster = new NEXeCluster(numXe, network, handlerRegistry);
	}
	assert(cluster != nullptr);

	// TODO when we have widespread C++14 support, use std::make_unique
	// and construct unique ptr and object pointed to in one memory operation.
	return std::unique_ptr<NECluster>(cluster);
}

std::unique_ptr<NECluster> NEClusterNetworkLoader::createNESuperCluster(
		int nTot, double numXe, double radius,
		IReactionNetwork& network) const {
	// Create the cluster
	auto superCluster = new NESuperCluster(numXe, nTot, nTot, radius, 0.0,
			network, handlerRegistry);

	// TODO when we have widespread C++14 support, use std::make_unique
	// and construct unique ptr and object pointed to in one memory operation.
	return std::unique_ptr<NECluster>(superCluster);
}

void NEClusterNetworkLoader::pushNECluster(
		std::unique_ptr<NEClusterReactionNetwork> & network,
		std::vector<std::reference_wrapper<Reactant> > & reactants,
		std::unique_ptr<NECluster> & cluster) {
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

std::unique_ptr<IReactionNetwork> NEClusterNetworkLoader::load(
		const IOptions& options) {
	// Get the dataset from the HDF5 files
	int normalSize = 0, superSize = 0;
	XFile networkFile(fileName);
	auto networkGroup = networkFile.getGroup<XFile::NetworkGroup>();
	assert(networkGroup);
	networkGroup->readNetworkSize(normalSize, superSize);

	// Initialization
	int numXe = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// Prepare the network
	std::unique_ptr<NEClusterReactionNetwork> network(
			new NEClusterReactionNetwork(handlerRegistry));

	// Loop on the clusters
	for (int i = 0; i < normalSize + superSize; i++) {
		// Open the cluster group
		XFile::ClusterGroup clusterGroup(*networkGroup, i);

		if (i < normalSize) {
			// Normal cluster
			// Read the composition
			auto comp = clusterGroup.readCluster(formationEnergy,
					migrationEnergy, diffusionFactor);
			numXe = comp[toCompIdx(Species::Xe)];

			// Create the cluster
			auto nextCluster = createNECluster(numXe, 0, 0, *network);

			// Set the formation energy
			nextCluster->setFormationEnergy(formationEnergy);
			// Set the diffusion factor and migration energy
			nextCluster->setMigrationEnergy(migrationEnergy);
			nextCluster->setDiffusionFactor(diffusionFactor);

			// Save it in the network
			pushNECluster(network, reactants, nextCluster);
		} else {
			// Super cluster
			int nTot = 0;
			double average = 0.0, radius = 0.0;
			clusterGroup.readNESuperCluster(nTot, average, radius);

			// Create the cluster
			auto nextCluster = createNESuperCluster(nTot, average, radius,
					*network);

			// Save it in the network
			pushNECluster(network, reactants, nextCluster);
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

std::unique_ptr<IReactionNetwork> NEClusterNetworkLoader::generate(
		const IOptions &options) {
	// Initial declarations
	int maxXe = options.getMaxImpurity();
	int numXe = 0;
	// Once we have C++14 support, use std::make_unique.
	std::unique_ptr<NEClusterReactionNetwork> network(
			new NEClusterReactionNetwork(handlerRegistry));
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// The diffusion factor for a single xenon in nm^2/s, 1.0 is used
	// to trigger the computation but the exact value and formula are
	// defined in NECluster
	double xeOneDiffusion = 1.0;
	// The migration energy for a single xenon in eV
	double xeOneMigration = 0.0;

	/**
	 * The set of xenon formation energies up to Xe_29 indexed by size. That is
	 * E_(f,Xe_1) = xeFormationEnergies[1]. The value at index zero is just
	 * padding to make the indexing easy.
	 */
	std::vector<double> xeFormationEnergies = { 7.0, 12.15, 17.15, 21.90, 26.50,
			31.05, 35.30, 39.45, 43.00, 46.90, 50.65, 53.90, 56.90, 59.80,
			62.55, 65.05, 67.45, 69.45, 71.20, 72.75, 74.15, 75.35, 76.40,
			77.25, 77.95, 78.45, 78.80, 78.95, 79.0 };

	// Generate the I clusters
	for (int i = 1; i <= maxXe; ++i) {
		// Set the composition
		numXe = i;
		// Create the cluster
		auto nextCluster = createNECluster(numXe, 0, 0, *network);

		// Set the other attributes
		if (i <= xeFormationEnergies.size())
			nextCluster->setFormationEnergy(xeFormationEnergies[i - 1]);
		else
			nextCluster->setFormationEnergy(79.0);
		if (i <= 1) {
			nextCluster->setDiffusionFactor(xeOneDiffusion);
			nextCluster->setMigrationEnergy(xeOneMigration);
		} else {
			nextCluster->setDiffusionFactor(0.0);
			nextCluster->setMigrationEnergy(
					std::numeric_limits<double>::infinity());
		}

		// Save it in the network
		pushNECluster(network, reactants, nextCluster);
	}

	// Set the network for all of the reactants. This MUST be done manually.
	for (Reactant& currCluster : reactants) {
		currCluster.updateFromNetwork();
	}

	// Create the reactions
	network->createReactionConnectivity();

	// Check if we want dummy reactions
	if (!dummyReactions) {
		// Apply sectional grouping
		applyGrouping(*network);
	}

//	// Dump the network we've created, if desired.
//	int rank;
//	auto xolotlComm = xolotlCore::MPIUtils::getMPIComm();
//	MPI_Comm_rank(xolotlComm, &rank);
//	if (rank == 0) {
//		// Dump the network we've created for comparison with baseline.
//		std::ofstream networkStream(netDebugOpts.second);
//		network->dumpTo(networkStream);
//	}

	// Need to use move() because return type uses smart pointer to base class,
	// not derived class that we created.
	// Some C++11 compilers accept it without the move, but apparently
	// that is not correct behavior until C++14.
	return std::move(network);
}

void NEClusterNetworkLoader::applyGrouping(IReactionNetwork& network) const {
	// Get the xenon cluster map
	auto const& xeMap = network.getAll(ReactantType::Xe);

	// Create a temporary vector for the loop
	std::vector<NECluster *> tempVector;

	// Initialize variables for the loop
	int count = 0, superCount = 0, width = sectionWidth;
	double size = 0.0, radius = 0.0, energy = 0.0;

	// Map to know which cluster is in which group
	std::map<int, int> clusterGroupMap;
	// Map to know which super cluster gathers which group
	// TODO replace use of naked ptr with ref.  Requires fix to use
	// refs in reactions instead of naked ptrs.
	std::map<int, NESuperCluster *> superGroupMap;

	// Loop on the xenon groups
	for (int k = xeMin; k <= xeMap.size(); k++) {
		// Get the corresponding cluster
		auto cluster = (NECluster *) network.get(Species::Xe, k);

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
		auto rawSuperCluster = new NESuperCluster(size, count, count, radius,
				energy, network, handlerRegistry);
		auto superCluster = std::unique_ptr<NESuperCluster>(rawSuperCluster);
		// Save access to the cluster so we can trigger updates
		// after we give it to the network.
		auto& scref = *superCluster;
		// Set the HeV vector
		scref.setXeVector(tempVector);
		// Give the cluster to the network.
		network.add(std::move(superCluster));

		// Keep the information of the group
		superGroupMap[superCount] = rawSuperCluster;

//		std::cout << superCount << " " << count << " "
//				<< rawSuperCluster->getName() << std::endl;

		// Reinitialize everything
		size = 0.0, radius = 0.0, energy = 0.0;
		count = 0;
		tempVector.clear();
		superCount++;
		width = max(sectionWidth * (int) (superCount / 10), sectionWidth);
	}

	// Tell each reactant to update the pairs vector with super clusters
	for (IReactant& currReactant : network.getAll()) {

		auto& cluster = static_cast<NECluster&>(currReactant);

		// Get their production and dissociation vectors
		// TODO can these be updated in place?
		auto react = cluster.reactingPairs;
		auto combi = cluster.combiningReactants;
		auto disso = cluster.dissociatingPairs;
		auto emi = cluster.emissionPairs;

		// Loop on its reacting pairs
		for (int l = 0; l < react.size(); l++) {
			// Test the first reactant
			if (react[l].first->getType() == ReactantType::Xe) {
				// Get its size
				auto nXe = react[l].first->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					auto newCluster = superGroupMap[clusterGroupMap[nXe]];
					react[l].first = newCluster;
					react[l].firstDistance = newCluster->getDistance(size);
				}
			}

			// Test the second reactant
			if (react[l].second->getType() == ReactantType::Xe) {
				// Get its size
				auto nXe = react[l].second->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					auto newCluster = superGroupMap[clusterGroupMap[nXe]];
					react[l].second = newCluster;
					react[l].secondDistance = newCluster->getDistance(nXe);
				}
			}
		}

		// Visit each combining reactant.
		std::for_each(combi.begin(), combi.end(),
				[this,&superGroupMap,&clusterGroupMap](NECluster::CombiningCluster& cc) {
					// Test the combining reactant
					NECluster& currCombining = *(cc.combining);
					if (currCombining.getType() == ReactantType::Xe) {
						// Get its size
						auto nXe = currCombining.getSize();
						// Test its size
						if (nXe >= xeMin) {
							// It has to be replaced by a super cluster
							auto newCluster = superGroupMap[clusterGroupMap[nXe]];
							cc.combining = newCluster;
							cc.distance = newCluster->getDistance(nXe);
						}
					}
				});

		// Loop on its dissociating pairs
		for (int l = 0; l < disso.size(); l++) {
			// Test the first reactant
			if (disso[l].first->getType() == ReactantType::Xe) {
				// Get its size
				auto nXe = disso[l].first->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					auto newCluster = superGroupMap[clusterGroupMap[nXe]];
					disso[l].first = newCluster;
					disso[l].firstDistance = newCluster->getDistance(nXe);
				}
			}

			// Test the second reactant
			if (disso[l].second->getType() == ReactantType::Xe) {
				// Get its size
				auto nXe = disso[l].second->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					auto newCluster = superGroupMap[clusterGroupMap[nXe]];
					disso[l].second = newCluster;
					disso[l].secondDistance = newCluster->getDistance(nXe);
				}
			}
		}

		// Loop on its emission pairs
		for (int l = 0; l < emi.size(); l++) {
			// Test the first reactant
			if (emi[l].first->getType() == ReactantType::Xe) {
				// Get its size
				auto nXe = emi[l].first->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					auto newCluster = superGroupMap[clusterGroupMap[nXe]];
					emi[l].first = newCluster;
					emi[l].firstDistance = newCluster->getDistance(nXe);
				}
			}

			// Test the second reactant
			if (emi[l].second->getType() == ReactantType::Xe) {
				// Get its size
				auto nXe = emi[l].second->getSize();
				// Test its size
				if (nXe >= xeMin) {
					// It has to be replaced by a super cluster
					auto newCluster = superGroupMap[clusterGroupMap[nXe]];
					emi[l].second = newCluster;
					emi[l].secondDistance = newCluster->getDistance(nXe);
				}
			}
		}

		// Set their production and dissociation vectors
		cluster.reactingPairs = react;
		cluster.combiningReactants = combi;
		cluster.dissociatingPairs = disso;
		cluster.emissionPairs = emi;
	}

	// Set the reaction network for each super reactant
	for (auto const& superMapItem : network.getAll(ReactantType::NESuper)) {
		auto& currCluster = static_cast<NESuperCluster&>(*(superMapItem.second));
		currCluster.updateFromNetwork();
	}

	// Remove Xe clusters bigger than xeMin from the network
	// Loop on the Xe clusters
	std::vector<std::reference_wrapper<IReactant> > doomedReactants;
	for (auto const& currMapItem : xeMap) {

		auto& currCluster = currMapItem.second;

		// Get the cluster's size.
		auto nXe = currCluster->getSize();

		// Check if the cluster is too large.
		if (nXe >= xeMin) {
			// The cluster is too large.  Add it to the ones we will remove.
			doomedReactants.push_back(*currCluster);
		}
	}
	network.removeReactants(doomedReactants);

	// Recompute Ids and network size
	network.reinitializeNetwork();

	return;
}

} // namespace xolotlCore

