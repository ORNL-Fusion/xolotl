// Includes
#include "HeInterstitialCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>

using namespace xolotlCore;

HeInterstitialCluster::HeInterstitialCluster(int numHelium, int numInterstitial) :
		PSICluster(1), numHe(numHelium), numI(numInterstitial) {
	size = numHe + numI;
	name = "HeI";
}

HeInterstitialCluster::HeInterstitialCluster(const HeInterstitialCluster &other) :
		PSICluster(other) {
	numHe = other.numHe;
	numI = other.numI;
}

HeInterstitialCluster::~HeInterstitialCluster() {
}

std::shared_ptr<Reactant> HeInterstitialCluster::clone() {
	std::shared_ptr<Reactant> reactant(new HeInterstitialCluster(*this));
	return reactant;
}

double HeInterstitialCluster::getGenByEm() {
	return 0;
}

double HeInterstitialCluster::getAnnByEm() {
	return 0;
}

int HeInterstitialCluster::getSpeciesSize(const std::string speciesName) {
	return 0;
}

void HeInterstitialCluster::createReactionConnectivity() {

	// Local Declarations
	auto psiNetwork = std::dynamic_pointer_cast<PSIClusterReactionNetwork>(
			network);
	auto props = psiNetwork->getProperties();
	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	int maxIClusterSize = std::stoi(props["maxIClusterSize"]);
	int maxHeVClusterSize = std::stoi(props["maxHeVClusterSize"]);
	int maxHeIClusterSize = std::stoi(props["maxHeIClusterSize"]);
	int index = 0;
	std::shared_ptr<Reactant> firstReactant, secondReactant;
	std::vector<int> firstComposition, secondComposition;

	// Connect this cluster to itself since any reaction will affect it
	index = network->getReactantId(*this) - 1;
	reactionConnectivity[index] = 1;

	/* ----- (He_a)(I_b) + (V_c) --> (He_a)[I_(b-c)] -----
	 * This section adds the clusters that produce this cluster to the array
	 * for vacancy absorption by HeI.
	 */
	for (int c = 1; c <= maxVClusterSize; c++) { // hehe... c++!
		// Get the first reactant's composition and then retrieve it
		firstComposition = psiNetwork->getCompositionVector(numHe, 0, numI + c);
		firstReactant = psiNetwork->getCompound("HeI", firstComposition);
		// Set the second reactant
		secondReactant = psiNetwork->get("V", c);
		// Create a ReactingPair with the two reactants
		if (firstReactant && secondReactant) {
			ReactingPair pair;
			pair.first = std::dynamic_pointer_cast<PSICluster>(firstReactant);
			pair.second = std::dynamic_pointer_cast<PSICluster>(secondReactant);
			// Add the pair to the list
			reactingPairs.push_back(pair);
		}
	}

	/* ----- (He_a)(I_b) + I --> (He_a)[I_(b+1)] -----
	 * This section adds the clusters that produce this cluster to the array
	 * for single species interstitial absorption by HeI.
	 *
	 * This section also handles the case of this cluster combining with a
	 * single interstitial to produce an HeI cluster with one addition
	 * interstitial.
	 */
	firstComposition = psiNetwork->getCompositionVector(numHe, 0, numI - 1);
	// Get those Reactants from the network
	firstReactant = psiNetwork->getCompound("HeI", firstComposition);
	secondReactant = psiNetwork->get("I", 1);
	if (firstReactant && secondReactant) {
		// Create the Reacting Pair
		ReactingPair pair;
		pair.first = std::dynamic_pointer_cast<PSICluster>(firstReactant);
		pair.second = std::dynamic_pointer_cast<PSICluster>(secondReactant);
		// Add the pair to the list
		reactingPairs.push_back(pair);
		// Add single I to the list of clusters this one interacts with if it
		// doesn't violate the maximum size limit.
		index = psiNetwork->getReactantId(*secondReactant) - 1;
		reactionConnectivity[index] = 1;
		combiningReactants.push_back(secondReactant);
	}

	/* ----- (He_a)(I_b) + He_c --> [He_(a+c)](I_b)
	 * HeI clusters can absorb helium clusters of any size so long as the
	 * maximum size limit is not violated.
	 */
	auto reactants = psiNetwork->getAll("He");
	combineClusters(reactants,maxHeIClusterSize,"HeI");

	/* ----- (A*He)(B*I) + (C*V) --> (A*He)[(B-C)*I] -----
	 * This section adds the clusters that are produced by this cluster to the
	 * array for vacancy absorption by HeI.
	 */
	reactants = psiNetwork->getAll("V");
	fillVWithI("V",reactants);

	return;
}

void HeInterstitialCluster::createDissociationConnectivity() {

	// Local Declarations
	auto psiNetwork = std::dynamic_pointer_cast<PSIClusterReactionNetwork>(
			network);
	auto props = psiNetwork->getProperties();
	int maxIClusterSize = std::stoi(props["maxIClusterSize"]), index = 0;
	std::vector<int> composition;
	std::shared_ptr<Reactant> singleCluster, otherMixedCluster;

	// He Dissociation
	composition = psiNetwork->getCompositionVector(numHe - 1, 0, numI);
	otherMixedCluster = psiNetwork->getCompound("HeI", composition);
	singleCluster = psiNetwork->get("He", 1);
	if (singleCluster && otherMixedCluster) {
		// Mark the reactant in the dissociation array
		index = psiNetwork->getReactantId(*otherMixedCluster) - 1;
		dissociationConnectivity[index] = 1;
		// Mark single He
		index = psiNetwork->getReactantId(*singleCluster) - 1;
		dissociationConnectivity[index] = 1;
	}

	// Interstitial Dissociation
	composition = psiNetwork->getCompositionVector(numHe, 0, numI - 1);
	otherMixedCluster = psiNetwork->getCompound("HeI", composition);
	singleCluster = psiNetwork->get("I", 1);
	if (singleCluster && otherMixedCluster) {
		// Mark the reactant in the dissociation array
		int index = psiNetwork->getReactantId(*otherMixedCluster) - 1;
		dissociationConnectivity[index] = 1;
		index = psiNetwork->getReactantId(*singleCluster) - 1;
		dissociationConnectivity[index] = 1;
	}

	return;
}

double HeInterstitialCluster::getDissociationFlux(double temperature) const {

	// Local Declarations
	std::map<std::string, int> composition;
	std::shared_ptr<PSICluster> currentReactant, secondReactant;
	double f4 = 0.0, f3 = 0.0;

	// Get the required dissociation clusters
	auto heCluster = std::dynamic_pointer_cast<PSICluster>(
			network->get("He", 1));
	auto vCluster = std::dynamic_pointer_cast<PSICluster>(network->get("V", 1));
	auto iCluster = std::dynamic_pointer_cast<PSICluster>(network->get("I", 1));

	// Only dissociate if possible
	if (heCluster && vCluster && iCluster) {
		// Calculate the much easier f4 term...
		f4 = calculateDissociationConstant(*this, *iCluster, temperature)
				+ calculateDissociationConstant(*this, *vCluster, temperature)
				+ calculateDissociationConstant(*this, *heCluster, temperature);

		// Loop over all the elements of the dissociation
		// connectivity to find where this mixed species dissociates
		auto reactants = network->getAll();
		int numReactants = dissociationConnectivity.size();
		for (int i = 0; i < numReactants; i++) {
			if (dissociationConnectivity[i] == 1) {
				// Set the current reactant
				currentReactant = std::dynamic_pointer_cast<PSICluster>(
						reactants->at(i));
				// Get the cluster map of this connection
				composition = currentReactant->getComposition();
				// We need to find if this is a Helium dissociation
				if (numHe - composition["He"] == 1 && numI == composition["I"]
						&& composition["V"] == 0) {
					secondReactant = heCluster;
				} else if (numHe == composition["He"]
						&& numI - composition["V"] == 1
						&& composition["V"] == 0) {
					// trap mutation
					secondReactant = vCluster;
				} else if (numHe == composition["He"]
						&& composition["I"] - numI == 1
						&& composition["V"] == 0) {
					// or interstitial dissociation
					secondReactant = iCluster;
				}
				// Update the flux calculation
				if (secondReactant) {
					f3 += calculateDissociationConstant(*currentReactant,
							*secondReactant, temperature)
							* currentReactant->getConcentration();
				}
			}
		}
	}

	return f3 - f4 * getConcentration();
}

std::map<std::string, int> HeInterstitialCluster::getClusterMap() {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = numHe;
	clusterMap["V"] = 0;
	clusterMap["I"] = numI;

	// Return it
	return clusterMap;
}

std::map<std::string, int> HeInterstitialCluster::getComposition() const {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = numHe;
	clusterMap["V"] = 0;
	clusterMap["I"] = numI;

	// Return it
	return clusterMap;
}

double HeInterstitialCluster::getReactionRadius() {
	return 0.0;
}
