/*
 * PSIClusterNetworkLoader.cpp
 *
 *  Created on: Mar 30, 2013
 *      Author: jaybilly
 */

#include <fstream>
#include "PSIClusterNetworkLoader.h"
#include <TokenizedLineReader.h>
#include <PSIHeCluster.h>
#include <PSIVCluster.h>
#include <PSIInterstitialCluster.h>
#include <PSIMixedCluster.h>
#include <PSIDCluster.h>
#include <PSITCluster.h>
#include <MathUtils.h>
#include <cassert>

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

std::unique_ptr<PSICluster> PSIClusterNetworkLoader::createPSICluster(int numHe,
		int numD, int numT, int numV, int numI,
		IReactionNetwork& network) const {

	// Local Declarations
	PSICluster* cluster = nullptr;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if ((numHe > 0 || numT > 0 || numD > 0) && numV > 0) {
		// Create a new MixedCluster
		cluster = new PSIMixedCluster(numHe, numD, numT, numV, network,
				handlerRegistry);
	} else if (numHe > 0 && numI > 0) {
		throw std::string("HeliumInterstitialCluster is not yet implemented.");
		// FIXME! Add code to add it to the list
	} else if (numHe > 0) {
		// Create a new HeCluster
		cluster = new PSIHeCluster(numHe, network, handlerRegistry);
	} else if (numV > 0) {
		// Create a new VCluster
		cluster = new PSIVCluster(numV, network, handlerRegistry);
	} else if (numI > 0) {
		// Create a new ICluster
		cluster = new PSIInterstitialCluster(numI, network, handlerRegistry);
	} else if (numD > 0) {
		// Create a new DCluster
		cluster = new PSIDCluster(numD, network, handlerRegistry);
	} else if (numT > 0) {
		// Create a new TCluster
		cluster = new PSITCluster(numT, network, handlerRegistry);
	}
	assert(cluster != nullptr);

	// TODO Once we have widespread C++14 support, use std::make_unique
	// instead of two steps (and two memory allocations).
	return std::unique_ptr<PSICluster>(cluster);
}

void PSIClusterNetworkLoader::pushPSICluster(
		std::unique_ptr<PSIClusterReactionNetwork> & network,
		std::vector<std::reference_wrapper<Reactant> > & reactants,
		std::unique_ptr<PSICluster> & cluster) {
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

PSIClusterNetworkLoader::PSIClusterNetworkLoader(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = nullptr;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;
	vMin = 1000000;
	sectionWidth[0] = 1;
	sectionWidth[1] = 1;
	sectionWidth[2] = 1;
	sectionWidth[3] = 1;
	maxHe = -1;
	maxI = -1;
	maxV = -1;
	maxD = -1;
	maxT = -1;

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
	sectionWidth[0] = 1;
	sectionWidth[1] = 1;
	sectionWidth[2] = 1;
	sectionWidth[3] = 1;
	maxHe = -1;
	maxI = -1;
	maxV = -1;
	maxD = -1;
	maxT = -1;

	return;
}

std::unique_ptr<IReactionNetwork> PSIClusterNetworkLoader::load(
		const IOptions& options) {
	// Local Declarations
	TokenizedLineReader<std::string> reader;
	std::vector<std::string> loadedLine;

	// TODO Once we have C++14, use std::make_unique.
	std::unique_ptr<PSIClusterReactionNetwork> network(
			new PSIClusterReactionNetwork(handlerRegistry));

	std::string error(
			"PSIClusterNetworkLoader Exception: Insufficient or erroneous data.");
	int numHe = 0, numV = 0, numI = 0, numW = 0, numD = 0, numT = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// Load the network if the stream is available
	if (networkStream != NULL) {
		// Load the stream
		reader.setInputStream(networkStream);

		// Loop over each line of the file, which should each be PSIClusters.
		loadedLine = reader.loadLine();
		while (loadedLine.size() > 0) {
			// Check the size of the loaded line
			if (loadedLine.size() != 8)
				// And notify the calling function if the size is insufficient
				throw error;
			// Load the sizes
			if (loadedLine[0][0] != '#') {
				numHe = std::stoi(loadedLine[0]);
				numD = std::stoi(loadedLine[1]);
				numT = std::stoi(loadedLine[2]);
				numV = std::stoi(loadedLine[3]);
				numI = std::stoi(loadedLine[4]);

				// If the cluster is big enough to be grouped
				if (numV >= vMin && numHe > 0) {
					// Created the coordinates and store them
					auto pair = std::make_tuple(numHe, numD, numT, numV);
					heVList.emplace(pair);
				} else {
					// Create the cluster
					auto nextCluster = createPSICluster(numHe, numD, numT, numV,
							numI, *network);
					// Load the energies
					formationEnergy = convertStrToDouble(loadedLine[5]);
					migrationEnergy = convertStrToDouble(loadedLine[6]);
					diffusionFactor = convertStrToDouble(loadedLine[7]);
					// Set the formation energy
					nextCluster->setFormationEnergy(formationEnergy);
					// Set the diffusion factor and migration energy
					nextCluster->setMigrationEnergy(migrationEnergy);
					nextCluster->setDiffusionFactor(diffusionFactor);

					// Save it in the network
					pushPSICluster(network, reactants, nextCluster);
				}
			}

			// Load the next line
			loadedLine = reader.loadLine();
		}

		// Update reactants now that they are in network.
		for (IReactant& currCluster : reactants) {
			currCluster.updateFromNetwork();
		}
	}

	// Check if we want dummy reactions
	if (!dummyReactions) {
		// Apply sectional grouping
		applySectionalGrouping(*network);
	}

	// Create the reactions
	network->createReactionConnectivity();

	// Recompute Ids and network size
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

std::unique_ptr<IReactionNetwork> PSIClusterNetworkLoader::generate(
		const IOptions &options) {
	// Initial declarations
	maxI = options.getMaxI(), maxHe = options.getMaxImpurity(), maxV =
			options.getMaxV(), maxD = options.getMaxD(), maxT =
			options.getMaxT();
	bool usePhaseCut = options.usePhaseCut();
	int numHe = 0, numD = 0, numT = 0, numV = 0, numI = 0;

	// Once we have C++14, use std::make_unique.
	std::unique_ptr<PSIClusterReactionNetwork> network(
			new PSIClusterReactionNetwork(handlerRegistry));
	std::vector<std::reference_wrapper<Reactant> > reactants;

	// I formation energies in eV
	std::vector<double> iFormationEnergies = { 10.0, 18.5, 27.0, 35.0, 42.5,
			48.0 };
	// I diffusion factors in nm^2/s
	std::vector<double> iDiffusion = { 8.8e+10, 8.0e+10, 3.9e+10, 2.0e+10,
			1.0e+10 };
	// I migration energies in eV
	std::vector<double> iMigration = { 0.01, 0.02, 0.03, 0.04, 0.05 };

	// He formation energies in eV
	std::vector<double> heFormationEnergies = { 6.15, 11.44, 16.35, 21.0, 26.1,
			30.24, 34.93, 38.80 };
	// He diffusion factors in nm^2/s
	std::vector<double> heDiffusion = { 2.9e+10, 3.2e+10, 2.3e+10, 1.7e+10,
			5.0e+09, 1.0e+09, 5.0e+08 };
	// He migration energies in eV
	std::vector<double> heMigration = { 0.13, 0.20, 0.25, 0.20, 0.12, 0.3, 0.4 };

	// The diffusion factor for a single deuterium.
	double dOneDiffusionFactor = 2.83e+11;
	// The migration energy for a single deuterium.
	double dOneMigrationEnergy = 0.38;

	// The diffusion factor for a single tritium.
	double tOneDiffusionFactor = 2.31e+11;
	// The migration energy for a single tritium.
	double tOneMigrationEnergy = 0.38;

	// The diffusion factor for a single vacancy in nm^2/s
	double vOneDiffusion = 1.8e+12;
	// The migration energy for a single vacancy in eV
	double vOneMigration = 1.30;

	/**
	 * The maximum number of helium atoms that can be combined with a vacancy
	 * cluster with size equal to the index i in the array plus one. For
	 * example, an HeV size cluster with size 1 would have size = i+1 = 1 and i
	 * = 0. It could support a mixture of up to nine helium atoms with one
	 * vacancy.
	 */
	std::vector<int> maxHePerV = { 9, 14, 18, 20, 27, 30, 35, 40, 45, 50, 55,
			60, 65, 70, 75, 80, 85, 90, 95, 98, 100, 101, 103, 105, 107, 109,
			110, 112, 116 };

	// Generate the I clusters
	for (int i = 1; i <= maxI; ++i) {
		// Set the composition
		numI = i;
		// Create the cluster
		auto nextCluster = createPSICluster(numHe, numD, numT, numV, numI,
				*network);

		// Set the other attributes
		if (i <= iFormationEnergies.size())
			nextCluster->setFormationEnergy(iFormationEnergies[i - 1]);
		else
			nextCluster->setFormationEnergy(48.0 + 6.0 * ((double) i - 6.0));
		if (i <= iDiffusion.size()) {
			nextCluster->setDiffusionFactor(iDiffusion[i - 1]);
			nextCluster->setMigrationEnergy(iMigration[i - 1]);
		} else {
			nextCluster->setDiffusionFactor(0.0);
			nextCluster->setMigrationEnergy(
					std::numeric_limits<double>::infinity());
		}

		// Save it in the network
		pushPSICluster(network, reactants, nextCluster);
	}

	// Reset the I composition
	numI = 0;

	// Generate the He clusters
	for (int i = 1; i <= maxHe; ++i) {
		// Set the composition
		numHe = i;
		// Create the cluster
		auto nextCluster = createPSICluster(numHe, numD, numT, numV, numI,
				*network);

		// Set the other attributes
		if (i <= heFormationEnergies.size())
			nextCluster->setFormationEnergy(heFormationEnergies[i - 1]);
		else
			nextCluster->setFormationEnergy(
					std::numeric_limits<double>::infinity());
		if (i <= heDiffusion.size()) {
			nextCluster->setDiffusionFactor(heDiffusion[i - 1]);
			nextCluster->setMigrationEnergy(heMigration[i - 1]);
		} else {
			nextCluster->setDiffusionFactor(0.0);
			nextCluster->setMigrationEnergy(
					std::numeric_limits<double>::infinity());
		}

		// Save it in the network
		pushPSICluster(network, reactants, nextCluster);
	}

	// Reset the He composition
	numHe = 0;

	// Generate the D clusters
	for (int i = 1; i <= maxD; ++i) {
		// Set the composition
		numD = i;
		// Create the cluster
		auto nextCluster = createPSICluster(numHe, numD, numT, numV, numI,
				*network);

		// Set the other attributes
		nextCluster->setFormationEnergy(
				std::numeric_limits<double>::infinity());
		if (i == 1) {
			nextCluster->setDiffusionFactor(dOneDiffusionFactor);
			nextCluster->setMigrationEnergy(dOneMigrationEnergy);
		} else {
			nextCluster->setDiffusionFactor(0.0);
			nextCluster->setMigrationEnergy(
					std::numeric_limits<double>::infinity());
		}

		// Save it in the network
		pushPSICluster(network, reactants, nextCluster);
	}

	// Reset the D composition
	numD = 0;

	// Generate the T clusters
	for (int i = 1; i <= maxT; ++i) {
		// Set the composition
		numT = i;
		// Create the cluster
		auto nextCluster = createPSICluster(numHe, numD, numT, numV, numI,
				*network);

		// Set the other attributes
		nextCluster->setFormationEnergy(
				std::numeric_limits<double>::infinity());
		if (i == 1) {
			nextCluster->setDiffusionFactor(tOneDiffusionFactor);
			nextCluster->setMigrationEnergy(tOneMigrationEnergy);
		} else {
			nextCluster->setDiffusionFactor(0.0);
			nextCluster->setMigrationEnergy(
					std::numeric_limits<double>::infinity());
		}

		// Save it in the network
		pushPSICluster(network, reactants, nextCluster);
	}

	// Reset the T composition
	numT = 0;

	// Check if the phase-cut method need to be applied
	if (usePhaseCut) {
		// Loop over vacancies in the outer loop.
		// This creates V and HeV.
		for (int i = 1; i <= maxV && i <= maxHePerV.size(); ++i) {
			// Create the V cluster
			numV = i;
			auto nextCluster = createPSICluster(numHe, numD, numT, numV, numI,
					*network);

			// Set its other attributes
			nextCluster->setFormationEnergy(getHeVFormationEnergy(numHe, numV));
			if (i <= 1) {
				nextCluster->setDiffusionFactor(vOneDiffusion);
				nextCluster->setMigrationEnergy(vOneMigration);
			} else {
				nextCluster->setDiffusionFactor(0.0);
				nextCluster->setMigrationEnergy(
						std::numeric_limits<double>::infinity());
			}

			// Save it in the network
			pushPSICluster(network, reactants, nextCluster);

			// For V < 12 loop on all the possible helium numbers up to
			// the maximum size in the maxHePerV array
			if (i < 12) {
				for (int j = 1; j <= maxHePerV[i - 1]; j++) {
					numHe = j;
					// Create the cluster
					nextCluster = createPSICluster(numHe, numD, numT, numV,
							numI, *network);
					// Set its attributes
					nextCluster->setFormationEnergy(
							getHeVFormationEnergy(numHe, numV));
					nextCluster->setDiffusionFactor(0.0);
					nextCluster->setMigrationEnergy(
							std::numeric_limits<double>::infinity());

					// Save it in the network
					pushPSICluster(network, reactants, nextCluster);
				}
			}
			// For bigger V only add the needed helium sizes
			else {
				for (int j = maxHePerV[i - 2] + 1; j <= maxHePerV[i - 1]; j++) {
					numHe = j;
					// Create the cluster
					nextCluster = createPSICluster(numHe, numD, numT, numV,
							numI, *network);
					// Set its attributes
					nextCluster->setFormationEnergy(
							getHeVFormationEnergy(numHe, numV));
					nextCluster->setDiffusionFactor(0.0);
					nextCluster->setMigrationEnergy(
							std::numeric_limits<double>::infinity());

					// Save it in the network
					pushPSICluster(network, reactants, nextCluster);
				}
			}

			// Reset the helium composition
			numHe = 0;
		}

		// Create V and HeV up to the maximum length with a constant nHe/nV
		// = 4,
		// Keeping only the last four ones
		for (int i = maxHePerV.size() + 1; i <= maxV; i++) {
			// Create the V cluster
			numV = i;
			auto nextCluster = createPSICluster(numHe, numD, numT, numV, numI,
					*network);

			// Set its other attributes
			nextCluster->setFormationEnergy(getHeVFormationEnergy(numHe, numV));
			if (i <= 1) {
				nextCluster->setDiffusionFactor(vOneDiffusion);
				nextCluster->setMigrationEnergy(vOneMigration);
			} else {
				nextCluster->setDiffusionFactor(0.0);
				nextCluster->setMigrationEnergy(
						std::numeric_limits<double>::infinity());
			}

			// Save it in the network
			pushPSICluster(network, reactants, nextCluster);

			// Loop on the helium number
			for (int j = (i * 4) - 3; j <= i * 4; j++) {
				numHe = j;
				// Create the cluster
				nextCluster = createPSICluster(numHe, numD, numT, numV, numI,
						*network);
				// Set its attributes
				nextCluster->setFormationEnergy(
						getHeVFormationEnergy(numHe, numV));
				nextCluster->setDiffusionFactor(0.0);
				nextCluster->setMigrationEnergy(
						std::numeric_limits<double>::infinity());

				// Save it in the network
				pushPSICluster(network, reactants, nextCluster);
			}

			// Reset the helium composition
			numHe = 0;
		}
	}
	// Else use the full network
	else {
		// Loop over vacancies in the outer loop.
		// This creates V and mixed clusters up to the maximum size in the
		// maxHePerV array.
		for (int i = 1; i <= maxV && i <= maxHePerV.size(); ++i) {
			// Create the V cluster
			numV = i;
			auto nextCluster = createPSICluster(numHe, numD, numT, numV, numI,
					*network);

			// Set its other attributes
			nextCluster->setFormationEnergy(getHeVFormationEnergy(numHe, numV));
			if (i <= 1) {
				nextCluster->setDiffusionFactor(vOneDiffusion);
				nextCluster->setMigrationEnergy(vOneMigration);
			} else {
				nextCluster->setDiffusionFactor(0.0);
				nextCluster->setMigrationEnergy(
						std::numeric_limits<double>::infinity());
			}

			// Save it in the network
			pushPSICluster(network, reactants, nextCluster);

			// Loop on the helium number
			int upperHe = maxHePerV[i - 1];
			if (maxHe <= 0)
				upperHe = 0;
			for (int j = 0; j <= upperHe; j++) {
				numHe = j;
				// Loop on the deuterium number
				int upperD = (int) ((2.0 / 3.0) * (double) maxHePerV[i - 1]);
				if (maxD <= 0)
					upperD = 0;
				for (int k = 0; k <= upperD; k++) {
					numD = k;
					// Loop on the tritium number
					int upperT = (int) ((2.0 / 3.0) * (double) maxHePerV[i - 1]);
					if (maxT <= 0)
						upperT = 0;
					for (int l = 0; l <= upperT; l++) {
						numT = l;

						if (numHe + numD + numT == 0)
							continue;

						// If the cluster is big enough to be grouped
						if (numV >= vMin) {
							// Created the coordinates and store them
							auto pair = std::make_tuple(numHe, numD, numT,
									numV);
							heVList.emplace(pair);
						} else {
							// Create the cluster
							nextCluster = createPSICluster(numHe, numD, numT,
									numV, numI, *network);
							// Set its attributes
							nextCluster->setFormationEnergy(
									getHeVFormationEnergy(numHe, numV));
							nextCluster->setDiffusionFactor(0.0);
							nextCluster->setMigrationEnergy(
									std::numeric_limits<double>::infinity());

							// Save it in the network
							pushPSICluster(network, reactants, nextCluster);
						}
					}
					// Reset the tritium composition
					numT = 0;
				}
				// Reset the deuterium composition
				numD = 0;
			}

			// Reset the helium composition
			numHe = 0;
		}

		// Create V and HeV up to the maximum length with a constant nHe/nV
		// = 4.
		for (int i = maxHePerV.size() + 1; i <= maxV; i++) {
			// Create the V cluster
			numV = i;
			auto nextCluster = createPSICluster(numHe, numD, numT, numV, numI,
					*network);

			// Set its other attributes
			nextCluster->setFormationEnergy(getHeVFormationEnergy(numHe, numV));
			if (i <= 1) {
				nextCluster->setDiffusionFactor(vOneDiffusion);
				nextCluster->setMigrationEnergy(vOneMigration);
			} else {
				nextCluster->setDiffusionFactor(0.0);
				nextCluster->setMigrationEnergy(
						std::numeric_limits<double>::infinity());
			}

			// Save it in the network
			pushPSICluster(network, reactants, nextCluster);

			// Loop on the helium number
			int upperHe = i * 4;
			if (maxHe <= 0)
				upperHe = 0;
			for (int j = 0; j <= upperHe; j++) {
				numHe = j;
				// Loop on the deuterium number
				int upperD = (int) ((2.0 / 3.0) * (double) i * 4.0);
				if (maxD <= 0)
					upperD = 0;
				for (int k = 0; k <= upperD; k++) {
					numD = k;
					// Loop on the tritium number
					int upperT = (int) ((2.0 / 3.0) * (double) i * 4.0);
					if (maxT <= 0)
						upperT = 0;
					for (int l = 0; l <= upperT; l++) {
						numT = l;

						if (numHe + numD + numT == 0)
							continue;

						// If the cluster is big enough to be grouped
						if (numV >= vMin) {
							// Created the coordinates and store them
							auto pair = std::make_tuple(numHe, numD, numT,
									numV);
							heVList.emplace(pair);
						} else {
							// Create the cluster
							nextCluster = createPSICluster(numHe, numD, numT,
									numV, numI, *network);
							// Set its attributes
							nextCluster->setFormationEnergy(
									getHeVFormationEnergy(numHe, numV));
							nextCluster->setDiffusionFactor(0.0);
							nextCluster->setMigrationEnergy(
									std::numeric_limits<double>::infinity());

							// Save it in the network
							pushPSICluster(network, reactants, nextCluster);
						}
					}
					// Reset the tritium composition
					numT = 0;
				}
				// Reset the deuterium composition
				numD = 0;
			}

			// Reset the helium composition
			numHe = 0;
		}
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

	// Recompute Ids and network size
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

double PSIClusterNetworkLoader::getHeVFormationEnergy(int numHe, int numV) {
	// V formation energies in eV
	std::array<double, 2> vFormationEnergies = { 3.6, 7.25 };

	// Coefficients for the Legendre polynomial fit
	// Low means V <= 27
	// Coefficients for c_0 in the 2D E_f,HeV fit
	std::array<double, 6> c0CoefficientsLow = { 253.35, 435.36, 336.50, 198.92,
			95.154, 21.544 };
	// Coefficients for c_1 in the 2D E_f,HeV fit
	std::array<double, 6> c1CoefficientsLow = { 493.29, 1061.3, 1023.9, 662.92,
			294.24, 66.962 };
	// Coefficients for c_2 in the 2D E_f,HeV fit
	std::array<double, 6> c2CoefficientsLow = { 410.40, 994.89, 1044.6, 689.41,
			286.52, 60.712 };
	// Coefficients for c_3 in the 2D E_f,HeV fit
	std::array<double, 6> c3CoefficientsLow = { 152.99, 353.16, 356.10, 225.75,
			87.077, 15.640 };
	// High means V > 27
	// Coefficients for c_0 in the 2D E_f,HeV fit
	std::array<double, 6> c0CoefficientsHigh = { -847.90, -3346.9, -4510.3,
			-3094.7, -971.18, -83.770 };
	// Coefficients for c_1 in the 2D E_f,HeV fit
	std::array<double, 6> c1CoefficientsHigh = { -1589.3, -4894.6, -6001.8,
			-4057.5, -1376.4, -161.91 };
	// Coefficients for c_2 in the 2D E_f,HeV fit
	std::array<double, 6> c2CoefficientsHigh = { 834.91, 1981.8, 1885.7, 1027.1,
			296.69, 29.902 };
	// Coefficients for c_3 in the 2D E_f,HeV fit
	std::array<double, 6> c3CoefficientsHigh = { 1547.2, 3532.3, 3383.6, 1969.2,
			695.17, 119.23 };

	/**
	 * The formation energies for He_xV_1. The entry at i = 0 is for a single
	 * vacancy (He_0V_1) and is there as a buffer. Like the formation energies,
	 * i = heSize.
	 */
	std::array<double, 14> heV1FormationEnergies = { 5.14166, 8.20919, 11.5304,
			14.8829, 18.6971, 22.2847, 26.3631, 30.1049, 34.0081, 38.2069,
			42.4217, 46.7378, 51.1551, 55.6738 };

	/**
	 * The formation energies for He_xV_2. The entry at i = 0 is for a
	 * di-vacancy (He_0V_2) and is there as a buffer. Like the formation
	 * energies, i = heSize.
	 */
	std::array<double, 18> heV2FormationEnergies =
			{ 7.10098, 8.39913, 9.41133, 11.8748, 14.8296, 17.7259, 20.7747,
					23.7993, 26.7984, 30.0626, 33.0385, 36.5173, 39.9406, 43.48,
					46.8537, 50.4484, 54.0879, 57.7939 };

	// Initial declarations
	double energy = -std::numeric_limits<double>::infinity();
	// The following coefficients are computed using the above and are used
	// to evaluate the full function f(x,y).
    std::array<double, 4> coefficients { 0.0, 0.0, 0.0, 0.0 };

	// Check to see if the vacancy size is large enough that the energy can
	// be computed from the fit or if it is so small that the exact values
	// must be used instead.
	if (numV > 2) {
		// Get the He/V ratio
		double x = 2.0 * (((double) numHe / (double) numV) / 9.0) - 1.0;
		// Initialize the vacancy number
		double y = 0.0;

		// We have 2 fits, one for low V and one for high V
		if (numV <= 27) {
			// Compute the vacancy number
			y = 2.0 * (((double) numV - 1.0) / 26.0) - 1.0;
			// Get the coefficients
			coefficients[0] = computeNthOrderLegendre<5>(x, c0CoefficientsLow);
			coefficients[1] = computeNthOrderLegendre<5>(x, c1CoefficientsLow);
			coefficients[2] = computeNthOrderLegendre<5>(x, c2CoefficientsLow);
			coefficients[3] = computeNthOrderLegendre<5>(x, c3CoefficientsLow);
		} else {
			// Compute the vacancy number
			y = 2.0 * (((double) numV - 1.0) / 451.0) - 1.0;
			// Get the coefficients
			coefficients[0] = computeNthOrderLegendre<5>(x, c0CoefficientsHigh);
			coefficients[1] = computeNthOrderLegendre<5>(x, c1CoefficientsHigh);
			coefficients[2] = computeNthOrderLegendre<5>(x, c2CoefficientsHigh);
			coefficients[3] = computeNthOrderLegendre<5>(x, c3CoefficientsHigh);
		}
		// Get the energy
		energy = computeNthOrderLegendre<3>(y, coefficients);

	} else if ((numV == 1 && numHe <= heV1FormationEnergies.size())
			|| (numV == 2 && numHe <= heV2FormationEnergies.size())) {
		// Get the exact energy
		energy =
				(numV == 1) ?
						heV1FormationEnergies[numHe - 1] :
						heV2FormationEnergies[numHe - 1];
	}
	// V Case
	if (numHe == 0 && numV < 3) {
		energy = (numV == 1) ? vFormationEnergies[0] : vFormationEnergies[1];

	}

	return energy;
}

void PSIClusterNetworkLoader::applySectionalGrouping(
		PSIClusterReactionNetwork& network) {

	// Define the phase space for the network
	int nDim = 1;
	Array<int, 5> list;
	list[0] = 0;
	// Add additional axis
	if (heVList.size() > 0) {
		if (maxHe > 0) {
			list[nDim] = 1;
			nDim++;
		}
		if (maxD > 0) {
			list[nDim] = 2;
			nDim++;
		}
		if (maxT > 0) {
			list[nDim] = 3;
			nDim++;
		}
		if (maxV > 0) {
			list[nDim] = 4;
			nDim++;
		}
	}

	// Initialize variables for the loop
	int count = 0, heIndex = ((maxHe > 0) && (maxD == 0) && (maxT == 0)),
			dIndex = ((maxD > 0) && (maxHe == 0) && (maxT == 0)), tIndex =
					((maxT > 0) && (maxHe == 0) && (maxD == 0)), vIndex = vMin,
			heWidth = sectionWidth[0], dWidth = sectionWidth[1], tWidth =
					sectionWidth[2], vWidth = sectionWidth[3],
			previousBiggestHe = 1, vMax = network.getMaxClusterSize(
					ReactantType::V), heMax = -1, dMax = -1, tMax = -1;
	double heSize = 0.0, dSize = 0.0, tSize = 0.0, vSize = 0.0;
	// Create a temporary vector for the loop to store which clusters go in which super cluster
	std::set<std::tuple<int, int, int, int> > tempVector;

	// Find the biggest helium, deuterium, tritium in the pairs
	for (auto pair : heVList) {
		if (std::get<0>(pair) > heMax)
			heMax = std::get<0>(pair);
		if (std::get<1>(pair) > dMax)
			dMax = std::get<1>(pair);
		if (std::get<2>(pair) > tMax)
			tMax = std::get<2>(pair);
	}

	// Take care of the clusters near the stability line
	// First get the biggest He cluster for the smaller V cluster
	IReactant::Composition comp;
	comp[toCompIdx(Species::He)] = previousBiggestHe;
	comp[toCompIdx(Species::V)] = vMin - 1;
	auto cluster = network.get(ReactantType::PSIMixed, comp);
	while (cluster) {
		previousBiggestHe++;
		comp[toCompIdx(Species::He)] = previousBiggestHe;
		cluster = network.get(ReactantType::PSIMixed, comp);
	}

	// Loop on the vacancy groups
	for (int k = vMin; k <= vMax; k++) {
		int heLow = previousBiggestHe, heHigh = -1, dLow = dMax, dHigh = -1,
				tLow = tMax, tHigh = -1;
		// Update the composition vector
		auto pair = std::make_tuple(previousBiggestHe, 0, 0, k);
		// While loop on the helium content because we don't know the upper bound
		while (heVList.find(pair) != heVList.end()) {
			// Will be used to know the actual widths of the group
			if (previousBiggestHe < heLow)
				heLow = previousBiggestHe;
			if (previousBiggestHe > heHigh)
				heHigh = previousBiggestHe;

			// Increment the counter
			count++;

			// Add this cluster to the temporary vector
			tempVector.emplace(pair);
			heSize += (double) previousBiggestHe;
			vSize += (double) k;
			// Remove the pair from the set because we don't need it anymore
			heVList.erase(pair);

			// Get the next cluster
			previousBiggestHe++;
			std::get<0>(pair) = previousBiggestHe;

			// Skip the last one
			if (previousBiggestHe == heMax)
				break;
		}

		// Check if there were clusters in this group
		if (count == 0)
			continue;

		// Average all values
		heSize = heSize / (double) count;
		vSize = vSize / (double) count;

		// Create the cluster
		double size[4] = { heSize, 0.0, 0.0, vSize };
		int width[4] = { heHigh - heLow + 1, 1, 1, 1 };
		PSISuperCluster* rawSuperCluster = new PSISuperCluster(size, count,
				width, network, handlerRegistry);
		auto superCluster = std::unique_ptr<PSISuperCluster>(rawSuperCluster);
		// Save access to the cluster so we can trigger updates
		// after we give it to the network.
		auto& scref = *superCluster;
		// Give the cluster to the network.
		network.add(std::move(superCluster));
		// Trigger cluster updates now it is in the network.
		scref.updateFromNetwork();
		// Set the HeV vector
		scref.setHeVVector(tempVector);

//		std::cout << "super: " << rawSuperCluster->getName() << " " << count
//				<< " 1" << std::endl;

		// Reinitialize everything
		heSize = 0.0, vSize = 0.0;
		count = 0;
		tempVector.clear();

		// Group the corresponding ones separately with D, T
		int nDGroup = dMax / sectionWidth[1] + 1;
		int nTGroup = tMax / sectionWidth[2] + 1;
		// Loop on the tritium groups
		for (int l = 0; l < nTGroup; l++) {
			// Loop on the deuterium groups
			for (int i = 0; i < nDGroup; i++) {
				// To compute the width
				int dLow = dMax, dHigh = -1, tLow = tMax, tHigh = -1;

				// Only the HeV clusters were grouped here, now group the corresponding HeDTV clusters
				// Loop on He, T, D
				for (int m = heLow; m <= heHigh; m++) {
					for (int o = tIndex; o < tIndex + tWidth; o++) {
						for (int n = dIndex; n < dIndex + dWidth; n++) {
							// Check if the corresponding coordinates are in the heVList set
							auto pair = std::make_tuple(m, n, o, k);
							if (heVList.find(pair) == heVList.end())
								continue;

							// Will be used to know the section width
							if (n < dLow)
								dLow = n;
							if (n > dHigh)
								dHigh = n;
							if (o < tLow)
								tLow = o;
							if (o > tHigh)
								tHigh = o;

							// Increment the counter
							count++;

							// Add this cluster coordinates to the temporary vector
							tempVector.emplace(pair);
							heSize += (double) m;
							dSize += (double) n;
							tSize += (double) o;
							// Remove the pair from the set because we don't need it anymore
							heVList.erase(pair);
						}
					}
				}

				// Check if there were clusters in this group
				if (count == 0) {
					// Reinitialize the group indices for the deuterium direction
					dIndex += dWidth;
					dWidth += sectionWidth[1];
					continue;
				}

				// Average all values
				heSize = heSize / (double) count;
				dSize = dSize / (double) count;
				tSize = tSize / (double) count;
				// Create the super cluster
				double size[4] = { heSize, dSize, tSize, (double) k };
				int width[4] = { heHigh - heLow + 1, dHigh - dLow + 1, tHigh
						- tLow + 1, 1 };
				PSISuperCluster* rawSuperCluster = new PSISuperCluster(size,
						count, width, network, handlerRegistry);

//				std::cout << "super: " << rawSuperCluster->getName() << " "
//						<< count << " " << width[0] << " " << width[1] << " "
//						<< width[2] << " " << width[3] << std::endl;

				auto superCluster = std::unique_ptr<PSISuperCluster>(
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
				heSize = 0.0, dSize = 0.0, tSize = 0.0;
				count = 0;
				tempVector.clear();
				// Reinitialize the group indices for the deuterium direction
				dIndex += dWidth;
				dWidth += sectionWidth[1];
			}

			// Reinitialize the group indices for the tritium direction
			tIndex += tWidth;
			dIndex = ((maxD > 0) && (maxHe == 0) && (maxT == 0));
			dWidth = sectionWidth[1];
		}

		tIndex = ((maxT > 0) && (maxHe == 0) && (maxD == 0));
	}

	// Get the number of groups in the helium and vacancy directions
	int nVGroup = (vMax - vMin) / sectionWidth[3] + 1;
	int nHeGroup = heMax / sectionWidth[0] + 1;
	int nDGroup = dMax / sectionWidth[1] + 1;
	int nTGroup = tMax / sectionWidth[2] + 1;

	// Loop on the vacancy groups
	for (int k = 0; k < nVGroup; k++) {
		// Loop on the tritium groups
		for (int l = 0; l < nTGroup; l++) {
			// Loop on the deuterium groups
			for (int i = 0; i < nDGroup; i++) {
				// Loop on the helium groups
				for (int j = 0; j < nHeGroup; j++) {
					// To compute the width
					int heLow = heMax, heHigh = -1, dLow = dMax,
							dHigh = -1, tLow = tMax, tHigh = -1, vLow = vMax,
							vHigh = -1;

					// Loop within the group
					for (int p = vIndex; p < vIndex + vWidth; p++) {
						for (int o = tIndex; o < tIndex + tWidth; o++) {
							for (int n = dIndex; n < dIndex + dWidth; n++) {
								for (int m = heIndex; m < heIndex + heWidth;
										m++) {
									// Check if the corresponding coordinates are in the heVList set
									auto pair = std::make_tuple(m, n, o, p);
									if (heVList.find(pair) == heVList.end())
										continue;

									// Skip the last ones
									if (m == heMax)
										continue;

									// Will be used to know the section width
									if (m < heLow)
										heLow = m;
									if (m > heHigh)
										heHigh = m;
									if (n < dLow)
										dLow = n;
									if (n > dHigh)
										dHigh = n;
									if (o < tLow)
										tLow = o;
									if (o > tHigh)
										tHigh = o;
									if (p < vLow)
										vLow = p;
									if (p > vHigh)
										vHigh = p;

									// Increment the counter
									count++;

									// Add this cluster coordinates to the temporary vector
									tempVector.emplace(pair);
									heSize += (double) m;
									dSize += (double) n;
									tSize += (double) o;
									vSize += (double) p;
									// Remove the pair from the set because we don't need it anymore
									heVList.erase(pair);
								}
							}
						}
					}

					// Check if there were clusters in this group
					if (count == 0) {
						// Reinitialize the group indices for the helium direction
						heIndex += heWidth;
						continue;
					}

					// Average all values
					heSize = heSize / (double) count;
					dSize = dSize / (double) count;
					tSize = tSize / (double) count;
					vSize = vSize / (double) count;
					// Create the super cluster
					double size[4] = { heSize, dSize, tSize, vSize };
					int width[4] = { heHigh - heLow + 1, dHigh - dLow + 1, tHigh
							- tLow + 1, vHigh - vLow + 1 };
					PSISuperCluster* rawSuperCluster = new PSISuperCluster(size,
							count, width, network, handlerRegistry);

//					std::cout << "big: " << rawSuperCluster->getName() << " "
//							<< count << " " << width[0] << " " << width[1]
//							<< " " << width[2] << " " << width[3] << std::endl;

					auto superCluster = std::unique_ptr<PSISuperCluster>(
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
					heSize = 0.0, dSize = 0.0, tSize = 0.0, vSize = 0.0;
					count = 0;
					tempVector.clear();
					// Reinitialize the group indices for the helium direction
					heIndex += heWidth;
				}

				// Reinitialize the group indices for the deuterium direction
				dIndex += dWidth;
				dWidth += sectionWidth[1];
				heIndex = (maxHe > 0) && (maxD == 0) && (maxT == 0);
			}

			// Reinitialize the group indices for the tritium direction
			tIndex += tWidth;
			dIndex = ((maxD > 0) && (maxHe == 0) && (maxT == 0));
			dWidth = sectionWidth[1];
		}

		// Reinitialize the group indices for the vacancy direction
		vIndex += vWidth;
		tIndex = ((maxT > 0) && (maxHe == 0) && (maxD == 0));
	}

	// Group the largest ones separately
	nDGroup = dMax / sectionWidth[1] + 1;
	nTGroup = tMax / sectionWidth[2] + 1;
	// Loop on the tritium groups
	for (int l = 0; l < nTGroup; l++) {
		// Loop on the deuterium groups
		for (int i = 0; i < nDGroup; i++) {
			// To compute the width
			int dLow = dMax, dHigh = -1, tLow = tMax, tHigh = -1;
			for (int o = tIndex; o < tIndex + tWidth; o++) {
				for (int n = dIndex; n < dIndex + dWidth; n++) {
					// Check if the corresponding coordinates are in the heVList set
					auto pair = std::make_tuple(heMax, n, o, vMax);
					if (heVList.find(pair) == heVList.end())
						continue;

					// Skip the last one
					if (o == tMax && n == dMax)
						continue;

					// Will be used to know the section width
					if (n < dLow)
						dLow = n;
					if (n > dHigh)
						dHigh = n;
					if (o < tLow)
						tLow = o;
					if (o > tHigh)
						tHigh = o;

					// Increment the counter
					count++;

					// Add this cluster coordinates to the temporary vector
					tempVector.emplace(pair);
					dSize += (double) n;
					tSize += (double) o;
					// Remove the pair from the set because we don't need it anymore
					heVList.erase(pair);
				}
			}

			// Check if there were clusters in this group
			if (count == 0) {
				// Reinitialize the group indices for the deuterium direction
				dIndex += dWidth;
				dWidth += sectionWidth[1];
				continue;
			}

			// Average all values
			dSize = dSize / (double) count;
			tSize = tSize / (double) count;
			// Create the super cluster
			double size[4] = { (double) heMax, dSize, tSize, (double) vMax };
			int width[4] = { 1, dHigh - dLow + 1, tHigh - tLow + 1, 1 };
			PSISuperCluster* rawSuperCluster = new PSISuperCluster(size, count,
					width, network, handlerRegistry);

//			std::cout << "last: " << rawSuperCluster->getName() << " " << count
//					<< " " << width[0] << " " << width[1] << " " << width[2]
//					<< " " << width[3] << std::endl;

			auto superCluster = std::unique_ptr<PSISuperCluster>(
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
			dSize = 0.0, tSize = 0.0;
			count = 0;
			tempVector.clear();
			// Reinitialize the group indices for the deuterium direction
			dIndex += dWidth;
			dWidth += sectionWidth[1];
		}

		// Reinitialize the group indices for the tritium direction
		tIndex += tWidth;
		dIndex = ((maxD > 0) && (maxHe == 0) && (maxT == 0));
	}

	// Do the very last one
	if (heVList.size() == 1) {
		auto pair = heVList.begin();
		// Create the super cluster corresponding to the biggest mixed ones
		// We know He and V already
		double size[4] = { (double) heMax, (double) dMax, (double) tMax, (double) vMax };
		tempVector.emplace(std::make_tuple(heMax, dMax, tMax, vMax));
		int width[4] = { 1, 1, 1, 1 };
		auto rawSuperCluster = new PSISuperCluster(size, 1, width, network,
				handlerRegistry);

		std::cout << "last: " << rawSuperCluster->getName() << std::endl;

		auto superCluster = std::unique_ptr<PSISuperCluster>(rawSuperCluster);
		// Save access to the cluster so we can trigger updates
		// after we give it to the network.
		auto& scref = *superCluster;
		// Give the cluster to the network.
		network.add(std::move(superCluster));
		// Trigger cluster updates now it is in the network.
		scref.updateFromNetwork();
		// Set the HeV vector
		scref.setHeVVector(tempVector);
	}

	// Now that all the clusters are created
	// Give the information on the phase space to the network
	network.setPhaseSpace(nDim, list);

	return;
}
