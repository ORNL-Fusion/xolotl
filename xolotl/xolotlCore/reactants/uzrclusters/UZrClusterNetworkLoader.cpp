#include "UZrClusterNetworkLoader.h"
#include <UZrXeCluster.h>
#include <UZrVCluster.h>
#include <UZrXeVCluster.h>
#include <UZrClusterReactionNetwork.h>
#include <xolotlPerf.h>
#include <MathUtils.h>
#include <cassert>
#include "xolotlCore/io/XFile.h"

namespace xolotlCore {

std::unique_ptr<UZrCluster> UZrClusterNetworkLoader::createUZrCluster(int numXe,
		int numV, IReactionNetwork &network) const {

	// Local Declarations
	UZrCluster *cluster = nullptr;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if (numXe > 0 && numV > 0) {
		// Create a new XeVCluster
		cluster = new UZrXeVCluster(numXe, numV, network, handlerRegistry);
	} else if (numXe > 0) {
		// Create a new XeCluster
		cluster = new UZrXeCluster(numXe, network, handlerRegistry);
	} else if (numV > 0) {
		// Create a new VCluster
		cluster = new UZrVCluster(numV, network, handlerRegistry);
	}
	assert(cluster != nullptr);

	return std::unique_ptr<UZrCluster>(cluster);
}

UZrClusterNetworkLoader::UZrClusterNetworkLoader(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = nullptr;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	maxXe = -1;
	maxV = -1;

	return;
}

UZrClusterNetworkLoader::UZrClusterNetworkLoader(
		const std::shared_ptr<std::istream> stream,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = nullptr;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	maxXe = -1;
	maxV = -1;

	return;
}

std::unique_ptr<IReactionNetwork> UZrClusterNetworkLoader::load(
		const IOptions &options) {
	// Get the dataset from the HDF5 files
	int normalSize = 0, superSize = 0;
	XFile networkFile(fileName);
	auto networkGroup = networkFile.getGroup<XFile::NetworkGroup>();
	assert(networkGroup);
	networkGroup->readNetworkSize(normalSize, superSize);

	// Initialization
	int numXe = 0, numV = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// Prepare the network
	std::unique_ptr<UZrClusterReactionNetwork> network(
			new UZrClusterReactionNetwork(handlerRegistry));

	// Set the lattice parameter in the network
	double latticeParam = options.getLatticeParameter();
	if (!(latticeParam > 0.0))
		latticeParam = uraniumZirconiumLatticeConstant;
	network->setLatticeParameter(latticeParam);

	// Set the xenon radius in the network
	double radius = options.getImpurityRadius();
	if (!(radius > 0.0))
		radius = xenonRadius;
	network->setImpurityRadius(radius);

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
			numV = comp[toCompIdx(Species::V)];

			// Create the cluster
			auto nextCluster = createUZrCluster(numXe, numV, *network);

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

			// Save access to it so we can trigger updates once
			// added to the network.
			reactants.emplace_back(*nextCluster);

			// Give the cluster to the network
			network->add(std::move(nextCluster));
		} else {
			// This is not happening yet
		}
	}

	// Ask reactants to update now that they are in network.
	for (IReactant &currReactant : reactants) {
		currReactant.updateFromNetwork();
	}

	// Set the reactions
	networkGroup->readReactions(*network);

	// Recompute Ids and network size
	network->reinitializeNetwork();

	return std::move(network);
}

std::unique_ptr<IReactionNetwork> UZrClusterNetworkLoader::generate(
		const IOptions &options) {
	// Initial declarations
	maxXe = options.getMaxImpurity(), maxV = options.getMaxV();
	int numXe = 0, numV = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;

	// Once we have C++14, use std::make_unique.
	std::unique_ptr<UZrClusterReactionNetwork> network(
			new UZrClusterReactionNetwork(handlerRegistry));
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// Set the lattice parameter in the network
	double latticeParam = options.getLatticeParameter();
	if (!(latticeParam > 0.0))
		latticeParam = uraniumZirconiumLatticeConstant;
	network->setLatticeParameter(latticeParam);

	// Set the xenon radius in the network
	double radius = options.getImpurityRadius();
	if (!(radius > 0.0))
		radius = xenonRadius;
	network->setImpurityRadius(radius);

	// TODO: replace with the correct kinetics here and phase-space

	// Xe formation energies in eV
	std::vector<double> heFormationEnergies = { 0.0 };
	// Xe diffusion factors in nm^2/s
	std::vector<double> heDiffusion = { 3.3e+10 };
	// Xe migration energies in eV
	std::vector<double> heMigration = { 0.06 };

	// V formation energies in eV
	std::vector<double> vFormationEnergies = { 0.0 };
	// V diffusion factors in nm^2/s
	std::vector<double> vDiffusion = { 2.5e+10 };
	// V migration energies in eV
	std::vector<double> vMigration = { 0.48 };

	// Generate the Xe clusters
	for (int i = 1; i <= maxXe; ++i) {
		// Set the composition
		numXe = i;
		// Create the cluster
		auto nextCluster = createUZrCluster(numXe, numV, *network);

		// Set the other attributes
		nextCluster->setFormationEnergy(0.0);
		if (i <= heDiffusion.size()) {
			nextCluster->setDiffusionFactor(heDiffusion[i - 1]);
			nextCluster->setMigrationEnergy(heMigration[i - 1]);
		} else {
			nextCluster->setDiffusionFactor(0.0);
			nextCluster->setMigrationEnergy(
					std::numeric_limits<double>::infinity());
		}

		// Save access to it so we can trigger updates once
		// added to the network.
		reactants.emplace_back(*nextCluster);

		// Give the cluster to the network
		network->add(std::move(nextCluster));
	}

	// Reset the Xe composition
	numXe = 0;

	// Generate the V clusters
	for (int i = 1; i <= maxV; ++i) {
		// Set the composition
		numV = i;
		// Create the cluster
		auto nextCluster = createUZrCluster(numXe, numV, *network);

		// Set the other attributes
		nextCluster->setFormationEnergy(0.0);
		if (i <= vDiffusion.size()) {
			nextCluster->setDiffusionFactor(vDiffusion[i - 1]);
			nextCluster->setMigrationEnergy(vMigration[i - 1]);
		} else {
			nextCluster->setDiffusionFactor(0.0);
			nextCluster->setMigrationEnergy(
					std::numeric_limits<double>::infinity());
		}

		// Save access to it so we can trigger updates once
		// added to the network.
		reactants.emplace_back(*nextCluster);

		// Give the cluster to the network
		network->add(std::move(nextCluster));
	}

	// Reset the V composition
	numV = 0;

	// Loop over vacancies in the outer loop.
	for (int i = 1; i <= maxV; ++i) {
		numV = i;

		// Loop on the xenon number
		for (int j = 1; j <= maxXe; j++) {
			numXe = j;

			// Create the cluster
			auto nextCluster = createUZrCluster(numXe, numV, *network);
			// Set its attributes
			nextCluster->setFormationEnergy(0.0);
			nextCluster->setDiffusionFactor(0.0);
			nextCluster->setMigrationEnergy(
					std::numeric_limits<double>::infinity());

			// Save access to it so we can trigger updates once
			// all are added to the network.
			reactants.emplace_back(*nextCluster);

			// Give the cluster to the network
			network->add(std::move(nextCluster));
		}
	}

	// Update reactants now that they are in network.
	for (IReactant &currCluster : reactants) {
		currCluster.updateFromNetwork();
	}

	// Create the reactions
	network->createReactionConnectivity();

	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeNetwork();

	return std::move(network);
}

} // namespace xolotlCore

