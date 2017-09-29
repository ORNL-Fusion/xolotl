/*
 * PSIClusterNetworkLoader.cpp
 *
 *  Created on: Mar 30, 2013
 *      Author: jaybilly
 */

#include "PSIClusterNetworkLoader.h"
#include <TokenizedLineReader.h>
#include <HeCluster.h>
#include <VCluster.h>
#include <InterstitialCluster.h>
#include <HeVCluster.h>
#include <PSISuperCluster.h>
// #include <HeInterstitialCluster.h>
#include <xolotlPerf.h>
#include <MathUtils.h>

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
		int numV, int numI) {
	// Local Declarations
	std::shared_ptr<PSICluster> cluster;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if (numHe > 0 && numV > 0) {
		// Create a new HeVCluster
		cluster = std::make_shared<HeVCluster>(numHe, numV, handlerRegistry);
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
	}

	return cluster;
}

void PSIClusterNetworkLoader::pushPSICluster(
		std::shared_ptr<PSIClusterReactionNetwork> & network,
		std::vector<std::shared_ptr<Reactant> > & reactants,
		std::shared_ptr<PSICluster> & cluster) {
	// Check if we want dummy reactions
	if (dummyReactions) {
		// Create a dummy cluster (Reactant) from the existing cluster
		auto dummyCluster = std::static_pointer_cast<Reactant>(
				cluster->Reactant::clone());
		// Add the cluster to the network
		network->add(dummyCluster);
		// Add it to the list so that we can set the network later
		reactants.push_back(dummyCluster);
	} else {
		// Add the cluster to the network
		network->add(cluster);
		// Add it to the list so that we can set the network later
		reactants.push_back(cluster);
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
	int numHe = 0, numV = 0, numI = 0;
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
			if (loadedLine.size() < 6)
				// And notify the calling function if the size is insufficient
				throw error;
			// Load the sizes
			if (loadedLine[0][0] != '#') {
				numHe = std::stoi(loadedLine[0]);
				numV = std::stoi(loadedLine[1]);
				numI = std::stoi(loadedLine[2]);
				// Create the cluster
				auto nextCluster = createPSICluster(numHe, numV, numI);
				// Load the energies
				formationEnergy = convertStrToDouble(loadedLine[3]);
				migrationEnergy = convertStrToDouble(loadedLine[4]);
				diffusionFactor = convertStrToDouble(loadedLine[5]);
				// Set the formation energy
				nextCluster->setFormationEnergy(formationEnergy);
				// Set the diffusion factor and migration energy
				nextCluster->setMigrationEnergy(migrationEnergy);
				nextCluster->setDiffusionFactor(diffusionFactor);

				// Push the cluster
				pushPSICluster(network, reactants, nextCluster);
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

std::shared_ptr<IReactionNetwork> PSIClusterNetworkLoader::generate(
		IOptions &options) {
	// Initial declarations
	int maxI = options.getMaxI(), maxHe = options.getMaxImpurity(), maxV =
			options.getMaxV();
	bool usePhaseCut = options.usePhaseCut();
	int numHe = 0, numV = 0, numI = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::shared_ptr<PSIClusterReactionNetwork> network = std::make_shared<
			PSIClusterReactionNetwork>(handlerRegistry);
	std::vector<std::shared_ptr<Reactant> > reactants;

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

	// V formation energies in eV
	std::vector<double> vFormationEnergies = { 3.6, 7.25 };
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
		auto nextCluster = createPSICluster(numHe, numV, numI);

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

		// Push the cluster
		pushPSICluster(network, reactants, nextCluster);
	}

	// Reset the I composition
	numI = 0;

	// Generate the He clusters
	for (int i = 1; i <= maxHe; ++i) {
		// Set the composition
		numHe = i;
		// Create the cluster
		auto nextCluster = createPSICluster(numHe, numV, numI);

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

		// Push the cluster
		pushPSICluster(network, reactants, nextCluster);
	}

	// Reset the He composition
	numHe = 0;

	// Check if the phase-cut method need to be applied
	if (usePhaseCut) {
		// Loop over vacancies in the outer loop.
		// This creates V and HeV.
		for (int i = 1; i <= maxV && i <= maxHePerV.size(); ++i) {
			// Create the V cluster
			numV = i;
			auto nextCluster = createPSICluster(numHe, numV, numI);

			// Set its other attributes
			if (i <= vFormationEnergies.size())
				nextCluster->setFormationEnergy(vFormationEnergies[i - 1]);
			else
				nextCluster->setFormationEnergy(getHeVFormationEnergy(0, i));
			if (i <= 1) {
				nextCluster->setDiffusionFactor(vOneDiffusion);
				nextCluster->setMigrationEnergy(vOneMigration);
			} else {
				nextCluster->setDiffusionFactor(0.0);
				nextCluster->setMigrationEnergy(
						std::numeric_limits<double>::infinity());
			}

			// Push the cluster
			pushPSICluster(network, reactants, nextCluster);

			// For V < 12 loop on all the possible helium numbers up to
			// the maximum size in the maxHePerV array
			if (i < 12) {
				for (int j = 1; j <= maxHePerV[i - 1]; j++) {
					numHe = j;
					// Create the cluster
					nextCluster = createPSICluster(numHe, numV, numI);
					// Set its attributes
					nextCluster->setFormationEnergy(
							getHeVFormationEnergy(j, i));
					nextCluster->setDiffusionFactor(0.0);
					nextCluster->setMigrationEnergy(
							std::numeric_limits<double>::infinity());

					// Push the cluster
					pushPSICluster(network, reactants, nextCluster);
				}
			}
			// For bigger V only add the needed helium sizes
			else {
				for (int j = maxHePerV[i - 2] + 1; j <= maxHePerV[i - 1]; j++) {
					numHe = j;
					// Create the cluster
					nextCluster = createPSICluster(numHe, numV, numI);
					// Set its attributes
					nextCluster->setFormationEnergy(
							getHeVFormationEnergy(j, i));
					nextCluster->setDiffusionFactor(0.0);
					nextCluster->setMigrationEnergy(
							std::numeric_limits<double>::infinity());

					// Push the cluster
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
			auto nextCluster = createPSICluster(numHe, numV, numI);

			// Set its other attributes
			if (i <= vFormationEnergies.size())
				nextCluster->setFormationEnergy(vFormationEnergies[i - 1]);
			else
				nextCluster->setFormationEnergy(getHeVFormationEnergy(0, i));
			if (i <= 1) {
				nextCluster->setDiffusionFactor(vOneDiffusion);
				nextCluster->setMigrationEnergy(vOneMigration);
			} else {
				nextCluster->setDiffusionFactor(0.0);
				nextCluster->setMigrationEnergy(
						std::numeric_limits<double>::infinity());
			}

			// Push the cluster
			pushPSICluster(network, reactants, nextCluster);

			// Loop on the helium number
			for (int j = (i * 4) - 3; j <= i * 4; j++) {
				numHe = j;
				// Create the cluster
				nextCluster = createPSICluster(numHe, numV, numI);
				// Set its attributes
				nextCluster->setFormationEnergy(getHeVFormationEnergy(j, i));
				nextCluster->setDiffusionFactor(0.0);
				nextCluster->setMigrationEnergy(
						std::numeric_limits<double>::infinity());

				// Push the cluster
				pushPSICluster(network, reactants, nextCluster);
			}

			// Reset the helium composition
			numHe = 0;
		}
	}
	// Else use the full network
	else {
		// Loop over vacancies in the outer loop.
		// This creates V and HeV up to the maximum size in the
		// maxHePerV array.
		for (int i = 1; i <= maxV && i <= maxHePerV.size(); ++i) {
			// Create the V cluster
			numV = i;
			auto nextCluster = createPSICluster(numHe, numV, numI);

			// Set its other attributes
			if (i <= vFormationEnergies.size())
				nextCluster->setFormationEnergy(vFormationEnergies[i - 1]);
			else
				nextCluster->setFormationEnergy(getHeVFormationEnergy(0, i));
			if (i <= 1) {
				nextCluster->setDiffusionFactor(vOneDiffusion);
				nextCluster->setMigrationEnergy(vOneMigration);
			} else {
				nextCluster->setDiffusionFactor(0.0);
				nextCluster->setMigrationEnergy(
						std::numeric_limits<double>::infinity());
			}

			// Push the cluster
			pushPSICluster(network, reactants, nextCluster);

			// Loop on the helium number
			for (int j = 1; j <= maxHePerV[i - 1]; j++) {
				numHe = j;
				// Create the cluster
				nextCluster = createPSICluster(numHe, numV, numI);
				// Set its attributes
				nextCluster->setFormationEnergy(getHeVFormationEnergy(j, i));
				nextCluster->setDiffusionFactor(0.0);
				nextCluster->setMigrationEnergy(
						std::numeric_limits<double>::infinity());

				// Push the cluster
				pushPSICluster(network, reactants, nextCluster);
			}

			// Reset the helium composition
			numHe = 0;
		}

		// Create V and HeV up to the maximum length with a constant nHe/nV
		// = 4.
		for (int i = maxHePerV.size() + 1; i <= maxV; i++) {
			// Create the V cluster
			numV = i;
			auto nextCluster = createPSICluster(numHe, numV, numI);

			// Set its other attributes
			if (i <= vFormationEnergies.size())
				nextCluster->setFormationEnergy(vFormationEnergies[i - 1]);
			else
				nextCluster->setFormationEnergy(getHeVFormationEnergy(0, i));
			if (i <= 1) {
				nextCluster->setDiffusionFactor(vOneDiffusion);
				nextCluster->setMigrationEnergy(vOneMigration);
			} else {
				nextCluster->setDiffusionFactor(0.0);
				nextCluster->setMigrationEnergy(
						std::numeric_limits<double>::infinity());
			}

			// Push the cluster
			pushPSICluster(network, reactants, nextCluster);

			// Loop on the helium number
			for (int j = 1; j <= i * 4; j++) {
				numHe = j;
				// Create the cluster
				nextCluster = createPSICluster(numHe, numV, numI);
				// Set its attributes
				nextCluster->setFormationEnergy(getHeVFormationEnergy(j, i));
				nextCluster->setDiffusionFactor(0.0);
				nextCluster->setMigrationEnergy(
						std::numeric_limits<double>::infinity());

				// Push the cluster
				pushPSICluster(network, reactants, nextCluster);
			}

			// Reset the helium composition
			numHe = 0;
		}
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
		applySectionalGrouping(network);
	}

	return network;
}

double PSIClusterNetworkLoader::getHeVFormationEnergy(int numHe, int numV) {
	// Coefficients for the Legendre polynomial fit
	// Low means V <= 27
	// Coefficients for c_0 in the 2D E_f,HeV fit
	std::vector<double> c0CoefficientsLow = { 253.35, 435.36, 336.50, 198.92,
			95.154, 21.544 };
	// Coefficients for c_1 in the 2D E_f,HeV fit
	std::vector<double> c1CoefficientsLow = { 493.29, 1061.3, 1023.9, 662.92,
			294.24, 66.962 };
	// Coefficients for c_2 in the 2D E_f,HeV fit
	std::vector<double> c2CoefficientsLow = { 410.40, 994.89, 1044.6, 689.41,
			286.52, 60.712 };
	// Coefficients for c_3 in the 2D E_f,HeV fit
	std::vector<double> c3CoefficientsLow = { 152.99, 353.16, 356.10, 225.75,
			87.077, 15.640 };
	// High means V > 27
	// Coefficients for c_0 in the 2D E_f,HeV fit
	std::vector<double> c0CoefficientsHigh = { -847.90, -3346.9, -4510.3,
			-3094.7, -971.18, -83.770 };
	// Coefficients for c_1 in the 2D E_f,HeV fit
	std::vector<double> c1CoefficientsHigh = { -1589.3, -4894.6, -6001.8,
			-4057.5, -1376.4, -161.91 };
	// Coefficients for c_2 in the 2D E_f,HeV fit
	std::vector<double> c2CoefficientsHigh = { 834.91, 1981.8, 1885.7, 1027.1,
			296.69, 29.902 };
	// Coefficients for c_3 in the 2D E_f,HeV fit
	std::vector<double> c3CoefficientsHigh = { 1547.2, 3532.3, 3383.6, 1969.2,
			695.17, 119.23 };

	/**
	 * The formation energies for He_xV_1. The entry at i = 0 is for a single
	 * vacancy (He_0V_1) and is there as a buffer. Like the formation energies,
	 * i = heSize.
	 */
	std::vector<double> heV1FormationEnergies = { 5.14166, 8.20919, 11.5304,
			14.8829, 18.6971, 22.2847, 26.3631, 30.1049, 34.0081, 38.2069,
			42.4217, 46.7378, 51.1551, 55.6738 };

	/**
	 * The formation energies for He_xV_2. The entry at i = 0 is for a
	 * di-vacancy (He_0V_2) and is there as a buffer. Like the formation
	 * energies, i = heSize.
	 */
	std::vector<double> heV2FormationEnergies =
			{ 7.10098, 8.39913, 9.41133, 11.8748, 14.8296, 17.7259, 20.7747,
					23.7993, 26.7984, 30.0626, 33.0385, 36.5173, 39.9406, 43.48,
					46.8537, 50.4484, 54.0879, 57.7939 };

	// Initial declarations
	double energy = -std::numeric_limits<double>::infinity();
	// The following coefficients are computed using the above and are used
	// to evaluate the full function f(x,y).
	std::vector<double> coefficients = { 0.0, 0.0, 0.0, 0.0 };

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
			coefficients[0] = compute5thOrderLegendre(x, c0CoefficientsLow);
			coefficients[1] = compute5thOrderLegendre(x, c1CoefficientsLow);
			coefficients[2] = compute5thOrderLegendre(x, c2CoefficientsLow);
			coefficients[3] = compute5thOrderLegendre(x, c3CoefficientsLow);
		} else {
			// Compute the vacancy number
			y = 2.0 * (((double) numV - 1.0) / 451.0) - 1.0;
			// Get the coefficients
			coefficients[0] = compute5thOrderLegendre(x, c0CoefficientsHigh);
			coefficients[1] = compute5thOrderLegendre(x, c1CoefficientsHigh);
			coefficients[2] = compute5thOrderLegendre(x, c2CoefficientsHigh);
			coefficients[3] = compute5thOrderLegendre(x, c3CoefficientsHigh);
		}
		// Get the energy
		energy = compute3rdOrderLegendre(y, coefficients);

	} else if ((numV == 1 && numHe <= heV1FormationEnergies.size())
			|| (numV == 2 && numHe <= heV2FormationEnergies.size())) {
		// Get the exact energy
		energy =
				(numV == 1) ?
						heV1FormationEnergies[numHe - 1] :
						heV2FormationEnergies[numHe - 1];
	}

	return energy;
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
	std::vector<int> compositionVector = { previousBiggestHe, vMin - 1, 0 };
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
		compositionVector[1] = k;
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
					std::vector<int> compositionVector = { m, n, 0 };
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
							composition[vType], 0 };
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
							composition[vType], 0 };
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
							composition[vType], 0 };
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
							composition[vType], 0 };
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
							composition[vType], 0 };
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
							composition[vType], 0 };
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
							composition[vType], 0 };
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
