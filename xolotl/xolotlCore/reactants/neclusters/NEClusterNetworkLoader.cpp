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
		int nTot, int maxXe, IReactionNetwork& network) const {
	// Create the cluster
	auto superCluster = new NESuperCluster(maxXe, nTot, network,
			handlerRegistry);

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
	xeMax = -1;
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
	xeMax = -1;
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

	// Set the lattice parameter in the network
	double latticeParam = options.getLatticeParameter();
	if (!(latticeParam > 0.0))
		latticeParam = uraniumDioxydeLatticeConstant;
	network->setLatticeParameter(latticeParam);

	// Set the xenon radius in the network
	double radius = options.getImpurityRadius();
	if (!(radius > 0.0))
		radius = xenonRadius;
	network->setImpurityRadius(radius);

	// Set the density in a bubble
	network->setDensity(options.getDensity());

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

			if (numXe == 1) {
				// If the diffusivity is given
				if (options.getXenonDiffusivity() > 0.0) {
					nextCluster->setDiffusionFactor(
							options.getXenonDiffusivity());
					nextCluster->setMigrationEnergy(-1.0);
				}
			}

			// Save it in the network
			pushNECluster(network, reactants, nextCluster);
		} else {
			// Super cluster
			int nTot = 0, maxXe = 0;
			clusterGroup.readNESuperCluster(nTot, maxXe);

			// Create the cluster
			auto nextCluster = createNESuperCluster(nTot, maxXe, *network);

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
	xeMax = options.getMaxImpurity();
	int numXe = 0;
	// Once we have C++14 support, use std::make_unique.
	std::unique_ptr<NEClusterReactionNetwork> network(
			new NEClusterReactionNetwork(handlerRegistry));
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// Set the lattice parameter in the network
	double latticeParam = options.getLatticeParameter();
	if (!(latticeParam > 0.0))
		latticeParam = uraniumDioxydeLatticeConstant;
	network->setLatticeParameter(latticeParam);

	// Set the xenon radius in the network
	double radius = options.getImpurityRadius();
	if (!(radius > 0.0))
		radius = xenonRadius;
	network->setImpurityRadius(radius);

	// Set the density in a bubble
	network->setDensity(options.getDensity());

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

	// Generate the Xe clusters
	for (int i = 1; i <= min(xeMin - 1, xeMax); ++i) {
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
			// If the diffusivity is given
			if (options.getXenonDiffusivity() > 0.0) {
				nextCluster->setDiffusionFactor(options.getXenonDiffusivity());
				nextCluster->setMigrationEnergy(-1.0);
			}
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

	// Check if we want dummy reactions
	if (!dummyReactions) {
		// Apply sectional grouping
		applyGrouping(*network);
	}

	// Create the reactions
	network->createReactionConnectivity();

	return std::move(network);
}

void NEClusterNetworkLoader::applyGrouping(IReactionNetwork& network) const {

	// Initialize variables for the loop
	int count = 0, superCount = 0, width = sectionWidth;
	int size = 0;

	// Loop on the xenon groups
	for (int k = xeMin; k < xeMax; k++) {

		// Increment the counter
		count++;

		// Track the size
		size = k;

		// Continue if we are not at the wanted width yet
		if (count < width && k < xeMax - 1)
			continue;

		// Create the cluster
		auto rawSuperCluster = new NESuperCluster(size, count, network,
				handlerRegistry);

//		std::cout << superCount << " " << count << " "
//				<< rawSuperCluster->getName() << std::endl;

		auto superCluster = std::unique_ptr<NESuperCluster>(rawSuperCluster);
		// Give the cluster to the network.
		network.add(std::move(superCluster));

		// Reinitialize everything
		size = 0;
		count = 0;
		superCount++;
		width += 1;
	}

	if (xeMin < xeMax) {
		// Group the last one alone
		auto rawSuperCluster = new NESuperCluster(xeMax, 1, network,
				handlerRegistry);

//		std::cout << superCount << " last " << rawSuperCluster->getName()
//				<< std::endl;

		auto superCluster = std::unique_ptr<NESuperCluster>(rawSuperCluster);
		// Give the cluster to the network.
		network.add(std::move(superCluster));
	}

	// Recompute Ids and network size
	network.reinitializeNetwork();

	return;
}

} // namespace xolotlCore

