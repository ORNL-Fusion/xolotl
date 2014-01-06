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

std::shared_ptr<PSICluster> HeInterstitialCluster::getThisSharedPtrFromNetwork() const {
	auto composition = getComposition();
	std::vector<int> compVec = { composition["He"], composition["V"],
			composition["I"] };
	return std::dynamic_pointer_cast<PSICluster>(
			network->getCompound(name, compVec));
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
	combineClusters(reactants, maxHeIClusterSize, "HeI");

	/* ----- (A*He)(B*I) + (C*V) --> (A*He)[(B-C)*I] -----
	 * This section adds the clusters that are produced by this cluster to the
	 * array for vacancy absorption by HeI.
	 */
	reactants = psiNetwork->getAll("V");
	fillVWithI("V", reactants);

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
	dissociateClusters(singleCluster, otherMixedCluster);

	// Interstitial Dissociation
	composition = psiNetwork->getCompositionVector(numHe, 0, numI - 1);
	otherMixedCluster = psiNetwork->getCompound("HeI", composition);
	singleCluster = psiNetwork->get("I", 1);
	dissociateClusters(singleCluster, otherMixedCluster);

	return;
}

double HeInterstitialCluster::getDissociationFlux(double temperature) const {

	// Local Declarations
	std::map<std::string, int> composition;
	std::shared_ptr<PSICluster> currentCluster, secondCluster;
	double f4 = 0.0, f3 = 0.0;

	// Get the required dissociating clusters
	auto heCluster = std::dynamic_pointer_cast<PSICluster>(
			network->get("He", 1));
	auto vCluster = std::dynamic_pointer_cast<PSICluster>(network->get("V", 1));
	auto iCluster = std::dynamic_pointer_cast<PSICluster>(network->get("I", 1));

	// Only dissociate if possible
	if (heCluster && vCluster && iCluster) {
		// FIXME! Make sure that this works as expected! Make sure that it
		// correctly picks out right component in
		// calculateDissociationConstant!
		// Calculate the much easier f4 term... first
		f4 = calculateDissociationConstant(*this, *heCluster, temperature)
				+ calculateDissociationConstant(*this, *vCluster, temperature)
				+ calculateDissociationConstant(*this, *iCluster, temperature);

		// Loop over all the elements of the dissociation
		// connectivity to find where this mixed species dissociates and
		// calculate the f3 term.
		auto reactants = network->getAll();
		int numClusters = dissociationConnectivity.size();
		for (int i = 0; i < numClusters; i++) {
			if (dissociationConnectivity[i] == 1) {
				// Set the current reactant
				currentCluster = std::dynamic_pointer_cast<PSICluster>(
						reactants->at(i));
				// Get the cluster map of this connection
				composition = currentCluster->getComposition();
				// We need to find if this is a Helium dissociation
				if (numHe - composition["He"] == 1 && numI == composition["I"]
						&& composition["V"] == 0) {
					secondCluster = heCluster;
				} else if (numHe == composition["He"]
						&& numI - composition["V"] == 1
						&& composition["V"] == 0) {
					// trap mutation
					secondCluster = vCluster;
				} else if (numHe == composition["He"]
						&& composition["I"] - numI == 1
						&& composition["V"] == 0) {
					// or interstitial dissociation
					secondCluster = iCluster;
				}
				// Update the flux calculation
				if (secondCluster) {
					// Get the single species cluster that comes out with it
					f3 += calculateDissociationConstant(*currentCluster,*secondCluster,
							temperature) * currentCluster->getConcentration();
				}
			}
		}
	}

	return f3 - (f4 * getConcentration());
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

double HeInterstitialCluster::getReactionRadius() const {
	return 0.0;
}

/**
 * This operation computes the partial derivatives due to dissociation
 * reactions. The partial derivatives due to dissociation for compound
 * clusters are significantly different than those single-species clusters.
 *
 * @param partials The vector into which the partial derivatives should be
 * inserted. This vector should have a length equal to the size of the
 * network.
 * @param temperature The temperature at which the reactions are occurring.
 */
void HeInterstitialCluster::getDissociationPartialDerivatives(
		std::vector<double> & partials, double temperature) const {

	// Get the required dissociation clusters
	auto heCluster = std::dynamic_pointer_cast<PSICluster>(
			network->get("He", 1));
	auto vCluster = std::dynamic_pointer_cast<PSICluster>(network->get("V", 1));
	auto iCluster = std::dynamic_pointer_cast<PSICluster>(network->get("I", 1));

	// Partial derivative with respect to changes in this cluster
	double partialDeriv = calculateDissociationConstant(*this, *heCluster,
			temperature)
			+ calculateDissociationConstant(*this, *vCluster, temperature)
			+ calculateDissociationConstant(*this, *iCluster, temperature);
	// Get the index
	int index = network->getReactantId(*getThisSharedPtrFromNetwork());
	// Add it to the list of partials
	partials[index] += partialDeriv;

	// Partial derivative with respect to an HeI cluster with one less helium
	std::vector<int> compositionVec = { numHe - 1, 0, numI };
	auto heIClusterLessHe = std::dynamic_pointer_cast<PSICluster>(
			network->getCompound("HeI", compositionVec));
	// Compute the partial derivative if the smaller clusters exists
	if (heIClusterLessHe) {
		partialDeriv = calculateDissociationConstant(*heIClusterLessHe,
				*heCluster, temperature);
		index = network->getReactantId(*heIClusterLessHe);
		partials[index] += partialDeriv;
	}

	// Partial derivative with respect to an HeI cluster with one less
	// vacancy
	compositionVec = { numHe, 0, numI - 1 };
	auto heIClusterLessI = std::dynamic_pointer_cast<PSICluster>(
			network->getCompound("HeI", compositionVec));
	// Compute the partial derivative if the smaller clusters exists
	if (heIClusterLessI) {
		partialDeriv = calculateDissociationConstant(*heIClusterLessI,
				*vCluster, temperature);
		index = network->getReactantId(*heIClusterLessI);
		partials[index] += partialDeriv;
	}

	// This cluster cannot dissociate into a smaller HeI cluster and a vacancy,
	// so there is no partial derivative term for that case.

	return;
}
