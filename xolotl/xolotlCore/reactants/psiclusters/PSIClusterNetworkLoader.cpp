/*
 * PSIClusterNetworkLoader.cpp
 *
 *  Created on: Mar 30, 2013
 *      Author: jaybilly
 */

#include "PSIClusterNetworkLoader.h"
#include <TokenizedLineReader.h>
#include <HeCluster.h>
#include <DCluster.h>
#include <TCluster.h>
#include <VCluster.h>
#include <InterstitialCluster.h>
#include <HeVCluster.h>
#include <PSISuperCluster.h>
// #include <HeInterstitialCluster.h>
#include <PSIClusterReactionNetwork.h>
#include <xolotlPerf.h>

using namespace xolotlCore;

/**
 * This operation converts a string to a double, taking in to account the fact
 * that the input file may contain keys such as "infinite."
 *
 * @param inString the string to be converted
 * @return the string as a double
 */
static inline double convertStrToDouble(const std::string& inString) {
	return (inString.compare("infinite") == 0) ?
			std::numeric_limits<double>::infinity() :
			strtod(inString.c_str(), NULL);
}

std::shared_ptr<PSICluster> PSIClusterNetworkLoader::createPSICluster(int numHe,
		int numV, int numI, int numD, int numT) {
	// Local Declarations
	std::shared_ptr<PSICluster> cluster;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if ((numHe > 0 || numT > 0 || numD > 0) && numV > 0) {
		// Create a new HeVCluster
		cluster = std::make_shared<HeVCluster>(numT, numD, numHe, numV, handlerRegistry);
	} else if (numHe > 0 && numI > 0) {
		throw std::string("HeliumInterstitialCluster is not yet implemented.");
		// FIXME! Add code to add it to the list
	} else if (numHe > 0) {
		// Create a new HeCluster
		cluster = std::make_shared<HeCluster>(numHe, handlerRegistry);
	} else if (numV > 0) {
		// Create a new VCluster
		cluster = std::make_shared<VCluster>(numV, handlerRegistry);
	} else if (numI > 0) {
		// Create a new ICluster
		cluster = std::make_shared<InterstitialCluster>(numI, handlerRegistry);
	} else if (numD > 0) {
		// Create a new DCluster
		cluster = std::make_shared<DCluster>(numD, handlerRegistry);
	} else if (numT > 0) {
		// Create a new TCluster
		cluster = std::make_shared<TCluster>(numT, handlerRegistry);
	}

	return cluster;
}

PSIClusterNetworkLoader::PSIClusterNetworkLoader(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = nullptr;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	vMin = 1000000;
	heSectionWidth = 1;
	vSectionWidth = 1;

	return;
}

PSIClusterNetworkLoader::PSIClusterNetworkLoader(
		const std::shared_ptr<std::istream> stream,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = nullptr;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	vMin = 1000000;
	heSectionWidth = 1;
	vSectionWidth = 1;

	return;
}

std::shared_ptr<IReactionNetwork> PSIClusterNetworkLoader::load() {
	// Local Declarations
	TokenizedLineReader<std::string> reader;
	std::vector<std::string> loadedLine;
	std::shared_ptr<PSIClusterReactionNetwork> network = std::make_shared<
			PSIClusterReactionNetwork>(handlerRegistry);

	std::string error(
			"PSIClusterNetworkLoader Exception: Insufficient or erroneous data.");
	int numHe = 0, numV = 0, numI = 0, numW = 0, numD = 0, numT = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::shared_ptr<Reactant> > reactants;

	// Load the network if the stream is available
	if (networkStream != NULL) {
		// Load the stream
		reader.setInputStream(networkStream);

		// Loop over each line of the file, which should each be PSIClusters.
		loadedLine = reader.loadLine();
		while (loadedLine.size() > 0) {
			// Check the size of the loaded line
			if (loadedLine.size() != 7)
				// And notify the calling function if the size is insufficient
				throw error;
			// Load the sizes
			if (loadedLine[0][0] != '#') {
				numHe = std::stoi(loadedLine[0]);
				numW = std::stoi(loadedLine[1]);
				numV = -numW * (numW < 0);
				numI = numW * (numW > 0);
				numD = std::stoi(loadedLine[2]);
				numT = std::stoi(loadedLine[3]);
				// Create the cluster
				auto nextCluster = createPSICluster(numHe, numV, numI, numD,
						numT);
				// Load the energies
				formationEnergy = convertStrToDouble(loadedLine[4]);
				migrationEnergy = convertStrToDouble(loadedLine[5]);
				diffusionFactor = convertStrToDouble(loadedLine[6]);
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

			// Load the next line
			loadedLine = reader.loadLine();
		}

		// Set the network for all of the reactants. This MUST be done manually.
		for (auto currCluster : reactants) {
			currCluster->setReactionNetwork(network);
		}
	}

	// Create the reactions
	network->createReactionConnectivity();

	// Check if we want dummy reactions
	if (!dummyReactions) {
		// Apply sectional grouping
		applySectionalGrouping(network);
	}

	return network;
}

void PSIClusterNetworkLoader::applySectionalGrouping(
		std::shared_ptr<IReactionNetwork> network) {
	// Get the HeV cluster map
	auto heVMap = network->getAll(heVType);

	// Create a temporary vector for the loop
	std::vector<PSICluster *> tempVector;

	// Initialize variables for the loop
	PSICluster * cluster = nullptr;
	std::shared_ptr<PSISuperCluster> superCluster;
	static std::map<std::string, int> composition;
	int count = 0, heIndex = 1, vIndex = vMin, heWidth = heSectionWidth,
			vWidth = vSectionWidth, previousBiggestHe = 1;
	double heSize = 0.0, vSize = 0.0, radius = 0.0, energy = 0.0;

	// Map to know which cluster is in which group
	std::map<std::vector<int>, std::pair<int, int> > clusterGroupMap;
	// Map to know which super cluster gathers which group
	std::map<std::pair<int, int>, PSISuperCluster *> superGroupMap;

	// Take care of the clusters near the stability line
	// First get the biggest He cluster for the smaller V cluster
	std::vector<int> compositionVector = { previousBiggestHe, -(vMin - 1) };
	cluster = (PSICluster *) network->getCompound(heVType, compositionVector);
	while (cluster) {
		previousBiggestHe++;
		compositionVector[0] = previousBiggestHe;
		cluster = (PSICluster *) network->getCompound(heVType,
				compositionVector);
	}
	// Loop on the vacancy groups
	for (int k = vMin; k <= network->getAll(vType).size(); k++) {
		// Update the composition vector
		compositionVector[0] = previousBiggestHe;
		compositionVector[1] = -k;
		cluster = (PSICluster *) network->getCompound(heVType,
				compositionVector);
		// While loop on the helium content because we don't know the upper bound
		while (cluster) {
			// Increment the counter
			count++;

			// Add this cluster to the temporary vector
			tempVector.push_back(cluster);
			heSize += (double) previousBiggestHe;
			vSize += (double) k;
			radius += cluster->getReactionRadius();
			energy += cluster->getFormationEnergy();
			// Keep the information of the group
			clusterGroupMap[compositionVector] = std::make_pair(-1, k);

			// Get the next cluster
			previousBiggestHe++;
			compositionVector[0] = previousBiggestHe;
			cluster = (PSICluster *) network->getCompound(heVType,
					compositionVector);
		}

		// Check if there were clusters in this group
		if (count == 0)
			continue;

		// Average all values
		heSize = heSize / (double) count;
		vSize = vSize / (double) count;
		radius = radius / (double) count;
		energy = energy / (double) count;
		// Create the cluster
		superCluster = std::make_shared<PSISuperCluster>(heSize, vSize, count,
				count, 1, radius, energy, handlerRegistry);
		// Set the HeV vector
		superCluster->setHeVVector(tempVector);
		// Add this cluster to the network and clusters
		network->addSuper(superCluster);
		// Keep the information of the group
		superGroupMap[std::make_pair(-1, k)] = superCluster.get();

//		std::cout << "super: " << superCluster->getName() << " " << count << " 1"
//				<< std::endl;

		// Reinitialize everything
		heSize = 0.0, vSize = 0.0, radius = 0.0, energy = 0.0;
		count = 0;
		tempVector.clear();
	}

	// Get the number of groups in the helium and vacancy directions
	int nVGroup = (network->getAll(vType).size() - vMin) / vSectionWidth + 1;
	int nHeGroup = previousBiggestHe / heSectionWidth + 1;

	// Loop on the vacancy groups
	for (int k = 0; k < nVGroup; k++) {
		// Loop on the helium groups
		for (int j = 0; j < nHeGroup; j++) {
			// To check if the group is full
			int heLow = previousBiggestHe, heHigh = -1, vLow = network->getAll(
					vType).size(), vHigh = -1;

			// Loop within the group
			for (int n = vIndex; n < vIndex + vWidth; n++) {
				for (int m = heIndex; m < heIndex + heWidth; m++) {
					// Get the corresponding cluster
					std::vector<int> compositionVector = { m, -n };
					// Get the product of the same type as the second reactant
					cluster = (PSICluster *) network->getCompound(heVType,
							compositionVector);

					// Verify if the cluster exists
					if (!cluster)
						continue;

					// Verify it was not already used
					if (clusterGroupMap.find(compositionVector)
							!= clusterGroupMap.end())
						continue;

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
					tempVector.push_back(cluster);
					heSize += (double) m;
					vSize += (double) n;
					radius += cluster->getReactionRadius();
					energy += cluster->getFormationEnergy();
					// Keep the information of the group
					clusterGroupMap[compositionVector] = std::make_pair(j, k);
				}
			}

			// Check if there were clusters in this group
			if (count == 0) {
				// Reinitialize the group indices for the helium direction
				heIndex += heWidth;
//				heWidth = std::max((int) std::pow((double) (j * heSectionWidth), 3.0) / 400000, heSectionWidth);
//				heWidth -= heWidth % heSectionWidth;
				continue;
			}

			// Average all values
			heSize = heSize / (double) count;
			vSize = vSize / (double) count;
			radius = radius / (double) count;
			energy = energy / (double) count;
			// Create the super cluster
			if (count == heWidth * vWidth) {
				// Everything is fine, the cluster is full
				superCluster = std::make_shared<PSISuperCluster>(heSize, vSize,
						count, heWidth, vWidth, radius, energy,
						handlerRegistry);

//				std::cout << "normal: " << superCluster->getName() << " "
//						<< heWidth << " " << vWidth
//						<< std::endl;
			} else {
				// The cluster is smaller than we thought because we are at the edge
				superCluster = std::make_shared<PSISuperCluster>(heSize, vSize,
						count, heHigh - heLow + 1, vHigh - vLow + 1, radius,
						energy, handlerRegistry);

//				std::cout << "irregular: " << superCluster->getName() << " "
//						<< heHigh - heLow + 1 << " " << vHigh - vLow + 1
//						<< std::endl;
			}
			// Set the HeV vector
			superCluster->setHeVVector(tempVector);
			// Add this cluster to the network and clusters
			network->addSuper(superCluster);
			// Keep the information of the group
			superGroupMap[std::make_pair(j, k)] = superCluster.get();

			// Reinitialize everything
			heSize = 0.0, vSize = 0.0, radius = 0.0, energy = 0.0;
			count = 0;
			tempVector.clear();
			// Reinitialize the group indices for the helium direction
			heIndex += heWidth;
//			heWidth = std::max((int) std::pow((double) (j * heSectionWidth), 3.0) / 400000, heSectionWidth);
//			heWidth -= heWidth % heSectionWidth;
		}

		// Reinitialize the group indices for the vacancy direction
		vIndex += vWidth;
//		vWidth = std::max((int) std::pow((double) (k * vSectionWidth), 3.0) / 100000, vSectionWidth);
//		vWidth -= vWidth % vSectionWidth;
//		heWidth = heSectionWidth;
		heIndex = 1;
	}

	// Initialize variables for the loop
	PSISuperCluster * newCluster;
	// Loop on all the reactants to update the pairs vector with super clusters
	auto reactants = network->getAll();
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
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType],
							-composition[vType] };
					newCluster =
							superGroupMap[clusterGroupMap[compositionVector]];
					react[l].first = newCluster;
					react[l].firstHeDistance = newCluster->getHeDistance(
							composition[heType]);
					react[l].firstVDistance = newCluster->getVDistance(
							composition[vType]);
				}
			}

			// Test the second reactant
			if (react[l].second->getType() == heVType) {
				// Get its composition
				composition = react[l].second->getComposition();
				// Test its size
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType],
							-composition[vType] };
					newCluster =
							superGroupMap[clusterGroupMap[compositionVector]];
					react[l].second = newCluster;
					react[l].secondHeDistance = newCluster->getHeDistance(
							composition[heType]);
					react[l].secondVDistance = newCluster->getVDistance(
							composition[vType]);
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
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType],
							-composition[vType] };
					newCluster =
							superGroupMap[clusterGroupMap[compositionVector]];
					combi[l].combining = newCluster;
					combi[l].heDistance = newCluster->getHeDistance(
							composition[heType]);
					combi[l].vDistance = newCluster->getVDistance(
							composition[vType]);
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
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType],
							-composition[vType] };
					newCluster =
							superGroupMap[clusterGroupMap[compositionVector]];
					disso[l].first = newCluster;
					disso[l].firstHeDistance = newCluster->getHeDistance(
							composition[heType]);
					disso[l].firstVDistance = newCluster->getVDistance(
							composition[vType]);
				}
			}

			// Test the second reactant
			if (disso[l].second->getType() == heVType) {
				// Get its composition
				composition = disso[l].second->getComposition();
				// Test its size
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType],
							-composition[vType] };
					newCluster =
							superGroupMap[clusterGroupMap[compositionVector]];
					disso[l].second = newCluster;
					disso[l].secondHeDistance = newCluster->getHeDistance(
							composition[heType]);
					disso[l].secondVDistance = newCluster->getVDistance(
							composition[vType]);
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
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType],
							-composition[vType] };
					newCluster =
							superGroupMap[clusterGroupMap[compositionVector]];
					emi[l].first = newCluster;
					emi[l].firstHeDistance = newCluster->getHeDistance(
							composition[heType]);
					emi[l].firstVDistance = newCluster->getVDistance(
							composition[vType]);
				}
			}

			// Test the second reactant
			if (emi[l].second->getType() == heVType) {
				// Get its composition
				composition = emi[l].second->getComposition();
				// Test its size
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType],
							-composition[vType] };
					newCluster =
							superGroupMap[clusterGroupMap[compositionVector]];
					emi[l].second = newCluster;
					emi[l].secondHeDistance = newCluster->getHeDistance(
							composition[heType]);
					emi[l].secondVDistance = newCluster->getVDistance(
							composition[vType]);
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
	auto superMap = network->getAll(PSISuperType);
	// Set the reaction network for each super reactant
	for (auto currCluster : superMap) {
		currCluster->setReactionNetwork(network);
	}

	// Remove HeV clusters bigger than vMin from the network
	// Loop on the HeV clusters
	std::vector<IReactant*> doomedReactants;
	for (auto currCluster : heVMap) {

		// Get the composition
		composition = currCluster->getComposition();

		// Check if the cluster is too large.
		if (composition[vType] >= vMin) {
			// The cluster is too large.  Add it to the ones we will remove.
			doomedReactants.push_back(currCluster);
		}
	}
	network->removeReactants(doomedReactants);

	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeNetwork();

	return;
}
