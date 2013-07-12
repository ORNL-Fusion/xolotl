// Includes
#include "VCluster.h"

using namespace xolotlCore;

VCluster::VCluster(int nV) :
		PSICluster(nV) {
	// Set the reactant name appropriately
	name = "Vacancy";
}

VCluster::~VCluster() {
}

std::vector<int> VCluster::getConnectivity() {

	// Local Declarations
	std::map<std::string, std::string> props = *(network->properties);
	int maxMixedSize = strtol(props["maxMixedClusterSize"].c_str(), NULL, 10);
	int numHe = strtol(props["numHeClusters"].c_str(), NULL, 10);
	int numV = strtol(props["numVClusters"].c_str(), NULL, 10);
	int numI = strtol(props["numIClusters"].c_str(), NULL, 10);
	int numMixed = strtol(props["numMixedClusters"].c_str(), NULL, 10);
	int clusterIndex = 0;
	int heSize = 1, mixedSize = 1;
	std::vector<int> connectivityArray(network->reactants->size(), 0);

	// Vacancies interact with everything except for vacancies bigger than they
	// would combine with to form vacancies larger than the size limit.

	// -----  A*He + B*V → (A*He)(B*V) -----
	// Vacancy clusters can interact with any helium cluster so long as the sum
	// of the number of helium atoms and vacancies does not produce a cluster
	// with a size greater than the maximum mixed-species cluster size.
	clusterIndex = 0;
	// Get the size of the first helium cluster
	heSize = (std::dynamic_pointer_cast<PSICluster>(
			network->reactants->at(clusterIndex)))->getSize();
	// Loop over all of the valid helium clusters
	while (size + heSize <= maxMixedSize && clusterIndex < numHe) {
		// Set the connectivity to 1
		connectivityArray.at(clusterIndex) = 1;
		// Increment the counter
		clusterIndex++;
		// Get the size of the next vacancy cluster
		heSize = (std::dynamic_pointer_cast<PSICluster>(
				network->reactants->at(clusterIndex)))->getSize();
	}

	//----- A*V + B*V --> (A+B)*V -----
	// This cluster should interact with all other clusters of the same type up
	// to the max size minus the size of this one to produce larger clusters.
	for (int i = numHe; i < numHe + numV - size; i++) {
		connectivityArray.at(i) = 1;
	}

	//----- A*I + B*V
	// → (A-B)*I, if A > B
	// → (B-I)*V, if A < B
	// → 0, if A = B -----
	// Vacancies always annihilate interstitials.
	for (int i = numHe + numV; i < numHe + numV + numI; i++) {
		connectivityArray.at(i) = 1;
	}

	// ----- (A*He)(B*V) + C*V → (A*He)[(B+C)*V] -----
	// Vacancies can interact with a mixed-species cluster so long as the sum of
	// the number of vacancy atoms and the size of the mixed-species cluster
	// does not exceed the maximum mixed-species cluster size.
	clusterIndex = numHe + numV + numI;
	// Get the size of the first vacancy cluster
	mixedSize = (std::dynamic_pointer_cast<PSICluster>(
			network->reactants->at(clusterIndex)))->getSize();
	// Loop over all of the mixed clusters -- FIX BOUNDS!
	while (size + mixedSize <= maxMixedSize
			&& clusterIndex < numHe + numV + numI + numMixed) {
		// Set the connectivity to 1
		connectivityArray.at(clusterIndex) = 1;
		// Increment the counter
		clusterIndex++;
		// Get the size of the next mixed cluster
		mixedSize = (std::dynamic_pointer_cast<PSICluster>(
				network->reactants->at(clusterIndex)))->getSize();
	}

	return connectivityArray;
}

double VCluster::getDissociationFlux(const double temperature) {

	// Local Declarations
	double diss = 0.0;
	int numV = 0, deltaIndex = -1;

	// Loop over all Reactants
	for (int j = 0; j < network->reactants->size(); j++) {
		numV = network->toClusterMap(j)["V"];
		// If the Jth reactant contains Vacancy, then we calculate
		if (numV > 0) {
			// Search for the index of the cluster that contains exactly
			// one less Vacancy than reactant->at(j)
			for (int k = 0; k < network->reactants->size(); k++) {
				if ((network->toClusterMap(k)["V"] - numV) == 1) {
					deltaIndex = k;
				}
			}

			// There may not have been an index that had one less
			// Vacancy, if so, we won't add to the dissociation flux
			if (deltaIndex != -1) {
				// Calculate the dissociation, with K^- evaluated
				// at deltaIndex and this Vacancy Cluster's index.
				diss = diss
						+ calculateDissociationConstant(j,
								network->toClusterIndex(getClusterMap()),
								temperature)
								* network->reactants->at(j)->getConcentration();
			}
		}
	}

	// Return the dissociation
	return diss;
}

bool VCluster::isProductReactant(int reactantI, int reactantJ) {

	// Local Declarations, integers for species number for I, J reactants
	int rI_I = 0, rJ_I = 0, rI_He = 0, rJ_He = 0, rI_V = 0, rJ_V = 0;

	// Get the ClusterMap corresponding to
	// the given reactants
	std::map<std::string, int> reactantIMap = network->toClusterMap(reactantI);
	std::map<std::string, int> reactantJMap = network->toClusterMap(reactantJ);

	// Grab the numbers for each species
	// from each Reactant
	rI_I = reactantIMap["I"];
	rJ_I = reactantJMap["I"];
	rI_He = reactantIMap["He"];
	rJ_He = reactantJMap["He"];
	rI_V = reactantIMap["V"];
	rJ_V = reactantJMap["V"];

	// We should have no interstitials, a
	// total of 0 Helium, and a total of
	// size Vacancies
	return ((rI_I + rJ_I) == 0) && ((rI_He + rJ_He) == 0)
			&& ((rI_V + rJ_V) == size);
}

std::map<std::string, int> VCluster::getClusterMap() {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = 0;
	clusterMap["V"] = size;
	clusterMap["I"] = 0;

	// Return it
	return clusterMap;
}
