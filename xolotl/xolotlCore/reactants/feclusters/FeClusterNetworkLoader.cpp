/*
 * FEClusterNetworkLoader.cpp
 */

#include "FeClusterNetworkLoader.h"
#include <FeHeCluster.h>
#include <FeVCluster.h>
#include <FeInterstitialCluster.h>
#include <FeHeVCluster.h>
#include <FeSuperCluster.h>
#include <FeClusterReactionNetwork.h>
#include <xolotlPerf.h>
#include <MathUtils.h>
#include <cassert>
#include "xolotlCore/io/XFile.h"

namespace xolotlCore {

std::unique_ptr<FeCluster> FeClusterNetworkLoader::createFeCluster(int numHe,
		int numV, int numI, IReactionNetwork& network) const {

	// Local Declarations
	FeCluster* cluster = nullptr;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if (numHe > 0 && numV > 0) {
		// Create a new HeVCluster
		cluster = new FeHeVCluster(numHe, numV, network, handlerRegistry);
	} else if (numHe > 0) {
		// Create a new HeCluster
		cluster = new FeHeCluster(numHe, network, handlerRegistry);
	} else if (numV > 0) {
		// Create a new VCluster
		cluster = new FeVCluster(numV, network, handlerRegistry);
	} else if (numI > 0) {
		// Create a new FeICluster
		cluster = new FeInterstitialCluster(numI, network, handlerRegistry);
	}
	assert(cluster != nullptr);

	// TODO Once we have widespread C++14 support, use std::make_unique
	// instead of two steps (and two memory allocations).
	return std::unique_ptr<FeCluster>(cluster);
}

FeClusterNetworkLoader::FeClusterNetworkLoader(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = nullptr;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	vMin = 1000000;
	heSectionWidth = 1;
	vSectionWidth = 1;
	maxHe = -1;
	maxI = -1;
	maxV = -1;

	return;
}

FeClusterNetworkLoader::FeClusterNetworkLoader(
		const std::shared_ptr<std::istream> stream,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = nullptr;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	vMin = 1000000;
	heSectionWidth = 1;
	vSectionWidth = 1;
	maxHe = -1;
	maxI = -1;
	maxV = -1;

	return;
}

std::unique_ptr<IReactionNetwork> FeClusterNetworkLoader::load(
		const IOptions& options) {

//	// Get the dataset from the HDF5 files
//	XFile networkFile(fileName);
//	auto networkGroup = networkFile.getGroup<XFile::NetworkGroup>();
//	assert(networkGroup);
//	auto networkVector = networkGroup->readNetwork();

	// Initialization
	int numHe = 0, numV = 0, numI = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// Prepare the network
	std::unique_ptr<FeClusterReactionNetwork> network(
			new FeClusterReactionNetwork(handlerRegistry));

//	// Loop on the networkVector
//	for (auto lineIt = networkVector.begin(); lineIt != networkVector.end();
//			lineIt++) {
//		// Composition of the cluster
//		numHe = (int) (*lineIt)[0];
//		numV = (int) (*lineIt)[1];
//		numI = (int) (*lineIt)[2];
//		// Create the cluster
//		auto nextCluster = createFeCluster(numHe, numV, numI, *network);
//
//		// Energies
//		formationEnergy = (*lineIt)[3];
//		migrationEnergy = (*lineIt)[4];
//		diffusionFactor = (*lineIt)[5];
//
//		// Set the formation energy
//		nextCluster->setFormationEnergy(formationEnergy);
//		// Set the diffusion factor and migration energy
//		nextCluster->setMigrationEnergy(migrationEnergy);
//		nextCluster->setDiffusionFactor(diffusionFactor);
//
//		// Check if we want dummy reactions
//		if (dummyReactions) {
//			// Create a dummy cluster (just a stock Reactant)
//			// from the existing cluster
//			// TODO Once C++11 support is widespread, use std::make_unique.
//			std::unique_ptr<Reactant> dummyCluster(new Reactant(*nextCluster));
//
//			// Keep a ref to it so we can trigger its updates after
//			// we add it to the network.
//			reactants.emplace_back(*dummyCluster);
//
//			// Give the cluster to the network
//			network->add(std::move(dummyCluster));
//
//		} else {
//			// Keep a ref to it so we can trigger its updates after
//			// we add it to the network.
//			reactants.emplace_back(*nextCluster);
//
//			// Give the cluster to the network
//			network->add(std::move(nextCluster));
//		}
//	}
//
//	// Ask reactants to update now that they are in network.
//	for (IReactant& currReactant : reactants) {
//		currReactant.updateFromNetwork();
//	}
//
//	// Check if we want dummy reactions
//	if (!dummyReactions) {
//		// Apply sectional grouping
//		applySectionalGrouping(*network);
//	}
//
//	// Create the reactions
//	network->createReactionConnectivity();
//
//	// Recompute Ids and network size and redefine the connectivities
//	network->reinitializeNetwork();

	// Need to use move() because return type uses smart pointer to base class,
	// not derived class that we created.
	// Some C++11 compilers accept it without the move, but apparently
	// that is not correct behavior until C++14.
	return std::move(network);
}

std::unique_ptr<IReactionNetwork> FeClusterNetworkLoader::generate(
		const IOptions &options) {
	// Initial declarations
	maxI = options.getMaxI(), maxHe = options.getMaxImpurity(), maxV =
			options.getMaxV();
	int numHe = 0, numV = 0, numI = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;

	// Once we have C++14, use std::make_unique.
	std::unique_ptr<FeClusterReactionNetwork> network(
			new FeClusterReactionNetwork(handlerRegistry));
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// I formation energies in eV
	std::vector<double> iFormationEnergies = { 0.0 };
	// I diffusion factors in nm^2/s
	std::vector<double> iDiffusion = { 1.0e+11 };
	// I migration energies in eV
	std::vector<double> iMigration = { 0.34 };

	// He formation energies in eV
	std::vector<double> heFormationEnergies = { 0.0 };
	// He diffusion factors in nm^2/s
	std::vector<double> heDiffusion = { 1.0e+11, 5.0e+10, 3.3e+10 };
	// He migration energies in eV
	std::vector<double> heMigration = { 0.06, 0.06, 0.06 };

	// V formation energies in eV
	std::vector<double> vFormationEnergies = { 0.0 };
	// V diffusion factors in nm^2/s
	std::vector<double> vDiffusion = { 1.0e+11, 5.0e+10, 3.3e+10, 2.5e+10 };
	// V migration energies in eV
	std::vector<double> vMigration = { 0.67, 0.62, 0.37, 0.48 };

	// Generate the I clusters
	for (int i = 1; i <= maxI; ++i) {
		// Set the composition
		numI = i;
		// Create the cluster
		auto nextCluster = createFeCluster(numHe, numV, numI, *network);

		// Set the other attributes
		nextCluster->setFormationEnergy(0.0);
		if (i <= iDiffusion.size()) {
			nextCluster->setDiffusionFactor(iDiffusion[i - 1]);
			nextCluster->setMigrationEnergy(iMigration[i - 1]);
		} else {
			nextCluster->setDiffusionFactor(0.0);
			nextCluster->setMigrationEnergy(
					std::numeric_limits<double>::infinity());
		}

		// Svae access to it so we can trigger updates once all are
		// added to the network.
		reactants.emplace_back(*nextCluster);

		// Give the cluster to the network
		network->add(std::move(nextCluster));
	}

	// Reset the I composition
	numI = 0;

	// Generate the He clusters
	for (int i = 1; i <= 8; ++i) {
		// Set the composition
		numHe = i;
		// Create the cluster
		auto nextCluster = createFeCluster(numHe, numV, numI, *network);

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

	// Reset the He composition
	numHe = 0;

	// Loop over vacancies in the outer loop.
	// This creates V and HeV up to the maximum size in the
	// maxHePerV array.
	for (int i = 1; i <= maxV; ++i) {
		// Create the V cluster
		numV = i;
		if (numV < 11) {
			auto nextCluster = createFeCluster(numHe, numV, numI, *network);

			// Set its other attributes
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
			// all are added to the network.
			reactants.emplace_back(*nextCluster);

			// Give the cluster to the network
			network->add(std::move(nextCluster));
		}

		// Loop on the helium number
		for (int j = 1; j <= maxHe; j++) {
			numHe = j;
			// Create the cluster only if it is not going to be grouped
			if (numHe < vMin && numV < vMin) {
				auto nextCluster = createFeCluster(numHe, numV, numI, *network);
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

		// Reset the helium composition
		numHe = 0;
	}

	// Update reactants now that they are in network.
	for (IReactant& currCluster : reactants) {
		currCluster.updateFromNetwork();
	}

	// Check if we want dummy reactions
	if (!dummyReactions) {
		// Apply sectional grouping
		applySectionalGrouping(*network);
	}

	// Create the reactions
	network->createReactionConnectivity();

	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeNetwork();

//	// Dump the network we've created, if desired.
//	int rank;
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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

void FeClusterNetworkLoader::applySectionalGrouping(
		FeClusterReactionNetwork& network) const {

	// Create a temporary vector for the loop
	std::vector<std::pair<int, int> > tempVector;

	// Initialize variables for the loop
	int count = 0, heIndex = 1, vIndex = 1, heWidth = heSectionWidth, vWidth =
			vSectionWidth;
	double heSize = 0.0, vSize = 0.0;

	// Get the number of groups in the helium and vacancy directions
	int nVGroup = maxV / vSectionWidth + 1;
	int nHeGroup = maxHe / heSectionWidth + 1;

	// Loop on the vacancy groups
	std::vector<IReactant::SizeType> superClusterBounds;
	for (int k = 0; k < nVGroup; k++) {
		// Add the bound the the network vector
		superClusterBounds.emplace_back(vIndex);

		// Loop on the helium groups
		for (int j = 0; j < nHeGroup; j++) {
			// To check if the group is full
			int heLow = maxHe, heHigh = -1, vLow = maxV, vHigh = -1;

			// Loop within the group
			for (int n = vIndex; n < vIndex + vWidth; n++) {
				if (n > maxV)
					continue;
				for (int m = heIndex; m < heIndex + heWidth; m++) {
					if (m > maxHe)
						continue;
					if (m < vMin && n < vMin)
						continue;
					// Get the corresponding cluster coordinates
					auto pair = std::make_pair(m, n);

					// Will be used to know if the group was full
					if (m < heLow)
						heLow = m;
					if (m > heHigh)
						heHigh = m;
					if (n < vLow)
						vLow = n;
					if (n > vHigh)
						vHigh = n;

					// Increment the counter
					count++;

					// Add this cluster to the temporary vector
					tempVector.push_back(pair);
					heSize += (double) m;
					vSize += (double) n;
				}
			}

			// Check if there were clusters in this group
			if (count == 0) {
				// Reinitialize the group indices for the helium direction
				heIndex += heWidth;
				heWidth = std::max(
						(int) std::pow((double) (j * heSectionWidth), 3.0)
								/ 4000, heSectionWidth);
				heWidth -= heWidth % heSectionWidth;
				continue;
			}

			// Average all values
			heSize = heSize / (double) count;
			vSize = vSize / (double) count;
			// Create the super cluster
			FeSuperCluster* rawSuperCluster = nullptr;
			if (count == heWidth * vWidth) {
				// Everything is fine, the cluster is full
				rawSuperCluster = new FeSuperCluster(heSize, vSize, count,
						heWidth, vWidth, network, handlerRegistry);

//				std::cout << "normal: " << rawSuperCluster->getName() << " "
//						<< heWidth << " " << vWidth << std::endl;
			} else {
				// The cluster is smaller than we thought because we are at the edge
				rawSuperCluster = new FeSuperCluster(heSize, vSize, count,
						heHigh - heLow + 1, vHigh - vLow + 1, network,
						handlerRegistry);

//				std::cout << "irregular: " << rawSuperCluster->getName() << " "
//						<< heHigh - heLow + 1 << " " << vHigh - vLow + 1
//						<< std::endl;
			}
			assert(rawSuperCluster != nullptr);
			auto superCluster = std::unique_ptr<FeSuperCluster>(
					rawSuperCluster);

			// Save access to the cluster so we can trigger updates
			// after we give it to the network.
			auto& scref = *superCluster;

			// Give the cluster to the network.
			network.add(std::move(superCluster));

			// Trigger cluster updates now it is in the network.
			scref.updateFromNetwork();

			// Set the HeV vector
			scref.setHeVVector(tempVector);

			// Reinitialize everything
			heSize = 0.0, vSize = 0.0;
			count = 0;
			tempVector.clear();
			// Reinitialize the group indices for the helium direction
			heIndex += heWidth;
			heWidth = std::max(
					(int) std::pow((double) (j * heSectionWidth), 3.0) / 4000,
					heSectionWidth);
			heWidth -= heWidth % heSectionWidth;

			if (heIndex > maxHe)
				break;
		}

		// Reinitialize the group indices for the vacancy direction
		vIndex += vWidth;
		vWidth = std::max(
				(int) std::pow((double) (k * vSectionWidth), 3.0) / 4000,
				vSectionWidth);
		vWidth -= vWidth % vSectionWidth;
		heWidth = heSectionWidth;
		heIndex = 1;

		if (vIndex > maxV)
			break;
	}

	// Add the bound the the network vector
	superClusterBounds.emplace_back(maxV + 1);

	// Now that we have the bound vector defined, tell the network to
	// build its lookup map for super clusters
	network.buildSuperClusterMap(superClusterBounds);

	return;
}

} // namespace xolotlCore

