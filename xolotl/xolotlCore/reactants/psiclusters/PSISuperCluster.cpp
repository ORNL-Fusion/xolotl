// Includes
#include "PSISuperCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

/**
 * The helium momentum partials.
 */
std::vector<double> heMomentumPartials;

/**
 * The vacancy momentum partials.
 */
std::vector<double> vMomentumPartials;

PSISuperCluster::PSISuperCluster(double numHe, double numV, int nTot,
		int heWidth, int vWidth, double radius, double energy,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(registry), numHe(numHe), numV(numV), nTot(nTot), l0(0.0), l1He(
				0.0), l1V(0.0), dispersionHe(0.0), dispersionV(0.0), heMomentumFlux(
				0.0), vMomentumFlux(0.0) {
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = (int) (numHe + numV);

	// Update the composition map
	compositionMap[heType] = (int) (numHe * (double) nTot);
	compositionMap[vType] = (int) (numV * (double) nTot);

	// Set the width
	sectionHeWidth = heWidth;
	sectionVWidth = vWidth;

	// Set the reaction radius and formation energy
	reactionRadius = radius;
	formationEnergy = 0.0; // It is set to 0.0 because we do not want the super clusters to undergo dissociation

	// Set the diffusion factor and the migration energy
	migrationEnergy = std::numeric_limits<double>::infinity();
	diffusionFactor = 0.0;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He_" << numHe << "V_" << numV;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = PSISuperType;

	return;
}

PSISuperCluster::PSISuperCluster(PSISuperCluster &other) :
		PSICluster(other) {
	numHe = other.numHe;
	numV = other.numV;
	nTot = other.nTot;
	sectionHeWidth = other.sectionHeWidth;
	sectionVWidth = other.sectionVWidth;
	l0 = other.l0;
	l1He = other.l1He;
	l1V = other.l1V;
	dispersionHe = other.dispersionHe;
	dispersionV = other.dispersionV;
	reactingMap = other.reactingMap;
	combiningMap = other.combiningMap;
	dissociatingMap = other.dissociatingMap;
	emissionMap = other.emissionMap;
	effReactingMap = other.effReactingMap;
	effCombiningMap = other.effCombiningMap;
	effDissociatingMap = other.effDissociatingMap;
	effEmissionMap = other.effEmissionMap;
	effReactingList = other.effReactingList;
	effCombiningList = other.effCombiningList;
	effDissociatingList = other.effDissociatingList;
	effEmissionList = other.effEmissionList;
	heMomentumFlux = other.heMomentumFlux;
	vMomentumFlux = other.vMomentumFlux;

	return;
}

std::shared_ptr<IReactant> PSISuperCluster::clone() {
	std::shared_ptr<IReactant> reactant(new PSISuperCluster(*this));

	return reactant;
}

double PSISuperCluster::getConcentration(double distHe, double distV) const {
	return l0 + (distHe * l1He) + (distV * l1V);
}

double PSISuperCluster::getHeMomentum() const {
	return l1He;
}

double PSISuperCluster::getVMomentum() const {
	return l1V;
}

double PSISuperCluster::getTotalConcentration() const {
	// Initial declarations
	int heIndex = 0, vIndex = 0;
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the vacancy width
	for (int k = 0; k < sectionVWidth; k++) {
		// Compute the vacancy index
		vIndex = (int) (numV - (double) sectionVWidth / 2.0) + k + 1;

		// Loop on the helium width
		for (int j = 0; j < sectionHeWidth; j++) {
			// Compute the helium index
			heIndex = (int) (numHe - (double) sectionHeWidth / 2.0) + j + 1;

			// Check if this cluster exists
			if (effReactingMap.find(std::make_pair(heIndex, vIndex))
					== effReactingMap.end())
				continue;

			// Compute the distances
			heDistance = getHeDistance(heIndex);
			vDistance = getVDistance(vIndex);

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, vDistance);
		}
	}

	return conc;
}

double PSISuperCluster::getTotalHeliumConcentration() const {
	// Initial declarations
	int heIndex = 0, vIndex = 0;
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the vacancy width
	for (int k = 0; k < sectionVWidth; k++) {
		// Compute the vacancy index
		vIndex = (int) (numV - (double) sectionVWidth / 2.0) + k + 1;

		// Loop on the helium width
		for (int j = 0; j < sectionHeWidth; j++) {
			// Compute the helium index
			heIndex = (int) (numHe - (double) sectionHeWidth / 2.0) + j + 1;

			// Check if this cluster exists
			if (effReactingMap.find(std::make_pair(heIndex, vIndex))
					== effReactingMap.end())
				continue;

			// Compute the distances
			heDistance = getHeDistance(heIndex);
			vDistance = getVDistance(vIndex);

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, vDistance) * (double) heIndex;
		}
	}

	return conc;
}

double PSISuperCluster::getTotalVacancyConcentration() const {
	// Initial declarations
	int heIndex = 0, vIndex = 0;
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the vacancy width
	for (int k = 0; k < sectionVWidth; k++) {
		// Compute the vacancy index
		vIndex = (int) (numV - (double) sectionVWidth / 2.0) + k + 1;

		// Loop on the helium width
		for (int j = 0; j < sectionHeWidth; j++) {
			// Compute the helium index
			heIndex = (int) (numHe - (double) sectionHeWidth / 2.0) + j + 1;

			// Check if this cluster exists
			if (effReactingMap.find(std::make_pair(heIndex, vIndex))
					== effReactingMap.end())
				continue;

			// Compute the distances
			heDistance = getHeDistance(heIndex);
			vDistance = getVDistance(vIndex);

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, vDistance) * (double) vIndex;
		}
	}

	return conc;
}

double PSISuperCluster::getHeDistance(int he) const {
	if (sectionHeWidth == 1)
		return 0.0;
	return 2.0 * (double) (he - numHe) / ((double) sectionHeWidth - 1.0);
}

double PSISuperCluster::getVDistance(int v) const {
	if (sectionVWidth == 1)
		return 0.0;
	return 2.0 * (double) (v - numV) / ((double) sectionVWidth - 1.0);
}

void PSISuperCluster::createReactionConnectivity() {
	// Aggregate the reacting pairs and combining reactants from the heVVector
	// Loop on the heVVector
	for (int i = 0; i < heVVector.size(); i++) {
		// Get the cluster composition
		auto comp = heVVector[i]->getComposition();
		// Get both production vectors
		auto react = heVVector[i]->reactingPairs;
		auto combi = heVVector[i]->combiningReactants;

		// Create the key to the map
		auto key = std::make_pair(comp[heType], comp[vType]);

		// Set them in the super cluster map
		reactingMap[key] = react;
		combiningMap[key] = combi;
	}

	return;
}

void PSISuperCluster::createDissociationConnectivity() {
	// Aggregate the dissociating and emission pairs from the heVVector
	// Loop on the heVVector
	for (int i = 0; i < heVVector.size(); i++) {
		// Get the cluster composition
		auto comp = heVVector[i]->getComposition();
		// Get both dissociation vectors
		auto disso = heVVector[i]->dissociatingPairs;
		auto emi = heVVector[i]->emissionPairs;

		// Create the key to the map
		auto key = std::make_pair(comp[heType], comp[vType]);

		// Set them in the super cluster map
		dissociatingMap[key] = disso;
		emissionMap[key] = emi;
	}

	return;
}

void PSISuperCluster::computeRateConstants() {
	// Local declarations
	PSICluster *firstReactant = nullptr, *secondReactant = nullptr,
			*combiningReactant = nullptr, *dissociatingCluster = nullptr,
			*otherEmittedCluster = nullptr, *firstCluster = nullptr,
			*secondCluster = nullptr;
	double rate = 0.0;
	int heIndex = 0, vIndex = 0;
	// Initialize the value for the biggest production rate
	double biggestProductionRate = 0.0;
	// Initialize the dispersion sum
	double nHeSquare = 0.0, nVSquare = 0.0;

	// Loop on the vacancy width
	for (int k = 0; k < sectionVWidth; k++) {
		// Compute the vacancy index
		vIndex = (int) (numV - (double) sectionVWidth / 2.0) + k + 1;

		// Loop on the helium width
		for (int j = 0; j < sectionHeWidth; j++) {
			// Compute the helium index
			heIndex = (int) (numHe - (double) sectionHeWidth / 2.0) + j + 1;

			// Create the key for the maps
			auto key = std::make_pair(heIndex, vIndex);

			// Check if this cluster exists
			if (reactingMap.find(key) == reactingMap.end())
				continue;

			// Compute nSquare for the dispersion
			nHeSquare += (double) heIndex * heIndex;
			nVSquare += (double) vIndex * vIndex;

			// Get all the reaction vectors at this index
			reactingPairs = reactingMap[key];
			combiningReactants = combiningMap[key];
			dissociatingPairs = dissociatingMap[key];
			emissionPairs = emissionMap[key];

			// Initialize all the effective vectors
			effReactingPairs.clear();
			effCombiningReactants.clear();
			effDissociatingPairs.clear();
			effEmissionPairs.clear();

			// Compute the reaction constant associated to the reacting pairs
			// Set the total number of reacting pairs
			int nPairs = reactingPairs.size();

			// Loop on them
			for (int i = 0; i < nPairs; i++) {
				// Get the reactants
				firstReactant = reactingPairs[i].first;
				secondReactant = reactingPairs[i].second;
				// Compute the reaction constant
				rate = calculateReactionRateConstant(*firstReactant,
						*secondReactant);
				// Set it in the pair
				reactingMap[key][i].kConstant = rate / (double) nTot;

				// Add the reacting pair to the effective vector
				// if the rate is not 0.0
				if (!xolotlCore::equal(rate, 0.0)) {
					effReactingPairs.push_back(&reactingMap[key][i]);

					// Check if the rate is the biggest one up to now
					if (rate > biggestProductionRate)
						biggestProductionRate = rate;
				}
			}

			// Compute the reaction constant associated to the combining reactants
			// Set the total number of combining reactants
			int nReactants = combiningReactants.size();
			// Loop on them
			for (int i = 0; i < nReactants; i++) {
				// Get the reactants
				combiningReactant = combiningReactants[i].combining;
				// Compute the reaction constant
				rate = calculateReactionRateConstant(*this, *combiningReactant);
				// Set it in the combining cluster
				combiningMap[key][i].kConstant = rate / (double) nTot;

				// Add the combining reactant to the effective vector
				// if the rate is not 0.0
				if (!xolotlCore::equal(rate, 0.0)) {
					effCombiningReactants.push_back(&combiningMap[key][i]);

					// Add itself to the list again to account for the correct rate
					if (id == combiningReactant->getId())
						effCombiningReactants.push_back(&combiningMap[key][i]);
				}
			}

			// Compute the dissociation constant associated to the dissociating clusters
			// Set the total number of dissociating clusters
			nPairs = dissociatingPairs.size();
			// Loop on them
			for (int i = 0; i < nPairs; i++) {
				dissociatingCluster = dissociatingPairs[i].first;
				// The second element of the pair is the cluster that is also
				// emitted by the dissociation
				otherEmittedCluster = dissociatingPairs[i].second;
				// Compute the dissociation constant
				// The order of the cluster is important here because of the binding
				// energy used in the computation. It is taken from the type of the first cluster
				// which must be the single one
				// otherEmittedCluster is the single size one
				rate = calculateDissociationConstant(*dissociatingCluster,
						*otherEmittedCluster, *this);

				// Set it in the pair
				dissociatingMap[key][i].kConstant = rate / (double) nTot;

				// Add the dissociating pair to the effective vector
				// if the rate is not 0.0
				if (!xolotlCore::equal(rate, 0.0)) {
					effDissociatingPairs.push_back(&dissociatingMap[key][i]);

					// Add itself to the list again to account for the correct rate
					if (id == otherEmittedCluster->getId())
						effDissociatingPairs.push_back(
								&dissociatingMap[key][i]);
				}
			}

			// Compute the dissociation constant associated to the emission of pairs of clusters
			// Set the total number of emission pairs
			nPairs = emissionPairs.size();
			// Loop on them
			for (int i = 0; i < nPairs; i++) {
				firstCluster = emissionPairs[i].first;
				secondCluster = emissionPairs[i].second;
				// Compute the dissociation rate
				rate = calculateDissociationConstant(*this, *firstCluster,
						*secondCluster);
				// Set it in the pair
				emissionMap[key][i].kConstant = rate / (double) nTot;

				// Add the emission pair to the effective vector
				// if the rate is not 0.0
				if (!xolotlCore::equal(rate, 0.0)) {
					effEmissionPairs.push_back(&emissionMap[key][i]);
				}
			}

			// Shrink the arrays to save some space
			effReactingPairs.shrink_to_fit();
			effCombiningReactants.shrink_to_fit();
			effDissociatingPairs.shrink_to_fit();
			effEmissionPairs.shrink_to_fit();

			// Set the arrays in the effective maps
			effReactingMap[key] = effReactingPairs;
			effCombiningMap[key] = effCombiningReactants;
			effDissociatingMap[key] = effDissociatingPairs;
			effEmissionMap[key] = effEmissionPairs;
		}
	}

	// Reset the vectors to save memory
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	// Set the biggest rate to the biggest production rate
	biggestRate = biggestProductionRate;

	// Compute the dispersions
	if (sectionHeWidth == 1)
		dispersionHe = 1.0;
	else
		dispersionHe = 2.0
				* (nHeSquare
						- ((double) compositionMap[heType]
								* ((double) compositionMap[heType]
										/ (double) nTot)))
				/ ((double) (nTot * (sectionHeWidth - 1)));

	if (sectionVWidth == 1)
		dispersionV = 1.0;
	else
		dispersionV = 2.0
				* (nVSquare
						- ((double) compositionMap[vType]
								* ((double) compositionMap[vType]
										/ (double) nTot)))
				/ ((double) (nTot * (sectionVWidth - 1)));

	// Method to optimize the reaction vectors
	optimizeReactions();

	return;
}

void PSISuperCluster::optimizeReactions() {
	// Local declarations
	double heFactor = 0.0, vFactor = 0.0, heDistance = 0.0, vDistance = 0.0;
	PSICluster *firstReactant = nullptr, *secondReactant = nullptr,
			*combiningReactant = nullptr, *dissociatingCluster = nullptr,
			*otherEmittedCluster = nullptr, *firstCluster = nullptr,
			*secondCluster = nullptr;
	int heIndex = 0, vIndex = 0;

	// Loop on the effective reacting map
	for (auto mapIt = effReactingMap.begin(); mapIt != effReactingMap.end();
			++mapIt) {
		// Get the pairs
		auto pairs = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = pairs.begin(); it != pairs.end();) {
			// Get the two reacting clusters
			firstReactant = (*it)->first;
			secondReactant = (*it)->second;

			// Create a new SuperClusterProductionPair
			SuperClusterProductionPair superPair(firstReactant, secondReactant,
					(*it)->kConstant);

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != effReactingMap.end();
					++mapItBis) {
				// Compute the helium index
				heIndex = mapItBis->first.first;
				heFactor = (double) (heIndex - numHe) / dispersionHe;
				// Compute the vacancy index
				vIndex = mapItBis->first.second;
				vFactor = (double) (vIndex - numV) / dispersionV;

				// Get the pairs
				auto pairsBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = pairsBis.begin(); itBis != pairsBis.end();) {
					// Get the two reacting clusters
					auto firstReactantBis = (*itBis)->first;
					auto secondReactantBis = (*itBis)->second;

					// Check if it is the same reaction
					if (firstReactantBis == firstReactant
							&& secondReactantBis == secondReactant) {
						// First is A, second is B, in A + B -> this
						superPair.a000 += 1.0;
						superPair.a001 += heFactor;
						superPair.a002 += vFactor;
						superPair.a100 += (*itBis)->firstHeDistance;
						superPair.a101 += (*itBis)->firstHeDistance * heFactor;
						superPair.a102 += (*itBis)->firstHeDistance * vFactor;
						superPair.a200 += (*itBis)->firstVDistance;
						superPair.a201 += (*itBis)->firstVDistance * heFactor;
						superPair.a202 += (*itBis)->firstVDistance * vFactor;
						superPair.a010 += (*itBis)->secondHeDistance;
						superPair.a011 += (*itBis)->secondHeDistance * heFactor;
						superPair.a012 += (*itBis)->secondHeDistance * vFactor;
						superPair.a020 += (*itBis)->secondVDistance;
						superPair.a021 += (*itBis)->secondVDistance * heFactor;
						superPair.a022 += (*itBis)->secondVDistance * vFactor;
						superPair.a110 += (*itBis)->firstHeDistance
								* (*itBis)->secondHeDistance;
						superPair.a111 += (*itBis)->firstHeDistance
								* (*itBis)->secondHeDistance * heFactor;
						superPair.a112 += (*itBis)->firstHeDistance
								* (*itBis)->secondHeDistance * vFactor;
						superPair.a120 += (*itBis)->firstHeDistance
								* (*itBis)->secondVDistance;
						superPair.a121 += (*itBis)->firstHeDistance
								* (*itBis)->secondVDistance * heFactor;
						superPair.a122 += (*itBis)->firstHeDistance
								* (*itBis)->secondVDistance * vFactor;
						superPair.a210 += (*itBis)->firstVDistance
								* (*itBis)->secondHeDistance;
						superPair.a211 += (*itBis)->firstVDistance
								* (*itBis)->secondHeDistance * heFactor;
						superPair.a212 += (*itBis)->firstVDistance
								* (*itBis)->secondHeDistance * vFactor;
						superPair.a220 += (*itBis)->firstVDistance
								* (*itBis)->secondVDistance;
						superPair.a221 += (*itBis)->firstVDistance
								* (*itBis)->secondVDistance * heFactor;
						superPair.a222 += (*itBis)->firstVDistance
								* (*itBis)->secondVDistance * vFactor;

						// Do not delete the element if it is the original one
						if (itBis == it) {
							++itBis;
							continue;
						}

						// Remove the reaction from the vector
						itBis = pairsBis.erase(itBis);
					}
					// Go to the next element
					else
						++itBis;
				}

				// Give back the pairs
				mapItBis->second = pairsBis;
			}

			// Add the super pair
			effReactingList.push_front(superPair);

			// Remove the reaction from the vector
			it = pairs.erase(it);
		}
	}

	// Loop on the effective combining map
	for (auto mapIt = effCombiningMap.begin(); mapIt != effCombiningMap.end();
			++mapIt) {
		// Get the pairs
		auto clusters = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = clusters.begin(); it != clusters.end();) {
			// Get the combining cluster
			combiningReactant = (*it)->combining;

			// Create a new SuperClusterProductionPair with NULL as the second cluster because
			// we do not need it
			SuperClusterProductionPair superPair(combiningReactant, NULL,
					(*it)->kConstant);

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != effCombiningMap.end();
					++mapItBis) {
				// Compute the helium index
				heIndex = mapItBis->first.first;
				heDistance = getHeDistance(heIndex);
				heFactor = (double) (heIndex - numHe) / dispersionHe;
				// Compute the vacancy index
				vIndex = mapItBis->first.second;
				vDistance = getVDistance(vIndex);
				vFactor = (double) (vIndex - numV) / dispersionV;

				// Get the pairs
				auto clustersBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = clustersBis.begin();
						itBis != clustersBis.end();) {
					// Get the two reacting clusters
					auto combiningReactantBis = (*itBis)->combining;

					// Check if it is the same reaction
					if (combiningReactantBis == combiningReactant) {
						// This is A, itBis is B, in this + B -> C
						superPair.a000 += 1.0;
						superPair.a001 += heFactor;
						superPair.a002 += vFactor;
						superPair.a010 += (*itBis)->heDistance;
						superPair.a011 += (*itBis)->heDistance * heFactor;
						superPair.a012 += (*itBis)->heDistance * vFactor;
						superPair.a020 += (*itBis)->vDistance;
						superPair.a021 += (*itBis)->vDistance * heFactor;
						superPair.a022 += (*itBis)->vDistance * vFactor;
						superPair.a100 += heDistance;
						superPair.a101 += heDistance * heFactor;
						superPair.a102 += heDistance * vFactor;
						superPair.a200 += vDistance;
						superPair.a201 += vDistance * heFactor;
						superPair.a202 += vDistance * vFactor;
						superPair.a110 += (*itBis)->heDistance * heDistance;
						superPair.a111 += (*itBis)->heDistance * heDistance
								* heFactor;
						superPair.a112 += (*itBis)->heDistance * heDistance
								* vFactor;
						superPair.a210 += (*itBis)->heDistance * vDistance;
						superPair.a211 += (*itBis)->heDistance * vDistance
								* heFactor;
						superPair.a212 += (*itBis)->heDistance * vDistance
								* vFactor;
						superPair.a120 += (*itBis)->vDistance * heDistance;
						superPair.a121 += (*itBis)->vDistance * heDistance
								* heFactor;
						superPair.a122 += (*itBis)->vDistance * heDistance
								* vFactor;
						superPair.a220 += (*itBis)->vDistance * vDistance;
						superPair.a221 += (*itBis)->vDistance * vDistance
								* heFactor;
						superPair.a222 += (*itBis)->vDistance * vDistance
								* vFactor;

						// Do not delete the element if it is the original one
						if (itBis == it) {
							++itBis;
							continue;
						}

						// Remove the reaction from the vector
						itBis = clustersBis.erase(itBis);
					}
					// Go to the next element
					else
						++itBis;
				}

				// Give back the pairs
				mapItBis->second = clustersBis;
			}

			// Add the super pair
			effCombiningList.push_front(superPair);

			// Remove the reaction from the vector
			it = clusters.erase(it);
		}
	}

	// Loop on the effective dissociating map
	for (auto mapIt = effDissociatingMap.begin();
			mapIt != effDissociatingMap.end(); ++mapIt) {
		// Get the pairs
		auto pairs = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = pairs.begin(); it != pairs.end();) {
			// Get the two reacting clusters
			dissociatingCluster = (*it)->first;
			otherEmittedCluster = (*it)->second;

			// Create a new SuperClusterProductionPair
			SuperClusterDissociationPair superPair(dissociatingCluster,
					otherEmittedCluster, (*it)->kConstant);

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != effDissociatingMap.end();
					++mapItBis) {
				// Compute the helium index
				heIndex = mapItBis->first.first;
				heFactor = (double) (heIndex - numHe) / dispersionHe;
				// Compute the vacancy index
				vIndex = mapItBis->first.second;
				vFactor = (double) (vIndex - numV) / dispersionV;

				// Get the pairs
				auto pairsBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = pairsBis.begin(); itBis != pairsBis.end();) {
					// Get the two reacting clusters
					auto dissociatingClusterBis = (*itBis)->first;
					auto otherEmittedClusterBis = (*itBis)->second;

					// Check if it is the same reaction
					if (dissociatingClusterBis == dissociatingCluster
							&& otherEmittedClusterBis == otherEmittedCluster) {
						// A is the dissociating cluster
						superPair.a00 += 1.0;
						superPair.a01 += heFactor;
						superPair.a02 += vFactor;
						superPair.a10 += (*itBis)->firstHeDistance;
						superPair.a11 += (*itBis)->firstHeDistance * heFactor;
						superPair.a12 += (*itBis)->firstHeDistance * vFactor;
						superPair.a20 += (*itBis)->firstVDistance;
						superPair.a21 += (*itBis)->firstVDistance * heFactor;
						superPair.a22 += (*itBis)->firstVDistance * vFactor;

						// Do not delete the element if it is the original one
						if (itBis == it) {
							++itBis;
							continue;
						}

						// Remove the reaction from the vector
						itBis = pairsBis.erase(itBis);
					}
					// Go to the next element
					else
						++itBis;
				}

				// Give back the pairs
				mapItBis->second = pairsBis;
			}

			// Add the super pair
			effDissociatingList.push_front(superPair);

			// Remove the reaction from the vector
			it = pairs.erase(it);
		}
	}

	// Loop on the effective emission map
	for (auto mapIt = effEmissionMap.begin(); mapIt != effEmissionMap.end();
			++mapIt) {
		// Get the pairs
		auto pairs = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = pairs.begin(); it != pairs.end();) {
			// Get the two reacting clusters
			firstCluster = (*it)->first;
			secondCluster = (*it)->second;

			// Create a new SuperClusterProductionPair
			SuperClusterDissociationPair superPair(firstCluster, secondCluster,
					(*it)->kConstant);

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != effEmissionMap.end();
					++mapItBis) {
				// Compute the helium index
				heIndex = mapItBis->first.first;
				heDistance = getHeDistance(heIndex);
				heFactor = (double) (heIndex - numHe) / dispersionHe;
				// Compute the vacancy index
				vIndex = mapItBis->first.second;
				vDistance = getVDistance(vIndex);
				vFactor = (double) (vIndex - numV) / dispersionV;

				// Get the pairs
				auto pairsBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = pairsBis.begin(); itBis != pairsBis.end();) {
					// Get the two reacting clusters
					auto firstClusterBis = (*itBis)->first;
					auto secondClusterBis = (*itBis)->second;

					// Check if it is the same reaction
					if (firstClusterBis == firstCluster
							&& secondClusterBis == secondCluster) {
						// A is the dissociating cluster
						superPair.a00 += 1.0;
						superPair.a01 += heFactor;
						superPair.a02 += vFactor;
						superPair.a10 += heDistance;
						superPair.a11 += heDistance * heFactor;
						superPair.a12 += heDistance * vFactor;
						superPair.a20 += vDistance;
						superPair.a21 += vDistance * heFactor;
						superPair.a22 += vDistance * vFactor;

						// Do not delete the element if it is the original one
						if (itBis == it) {
							++itBis;
							continue;
						}

						// Remove the reaction from the vector
						itBis = pairsBis.erase(itBis);
					}
					// Go to the next element
					else
						++itBis;
				}

				// Give back the pairs
				mapItBis->second = pairsBis;
			}

			// Add the super pair
			effEmissionList.push_front(superPair);

			// Remove the reaction from the vector
			it = pairs.erase(it);
		}
	}

	// Clear the maps because they won't be used anymore
	effReactingPairs.clear();
	effCombiningReactants.clear();
	effDissociatingPairs.clear();
	effEmissionPairs.clear();
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	return;
}

void PSISuperCluster::updateRateConstants() {
	// Local declarations
	PSICluster *firstReactant = nullptr, *secondReactant = nullptr,
			*combiningReactant = nullptr, *dissociatingCluster = nullptr,
			*otherEmittedCluster = nullptr, *firstCluster = nullptr,
			*secondCluster = nullptr;
	double rate = 0.0;
	// Initialize the value for the biggest production rate
	double biggestProductionRate = 0.0;

	// Loop on the reacting list
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
		// Get the reactants
		firstReactant = (*it).first;
		secondReactant = (*it).second;
		// Compute the reaction constant
		rate = calculateReactionRateConstant(*firstReactant, *secondReactant);
		// Set it in the pair
		(*it).kConstant = rate / (double) nTot;

		// Check if the rate is the biggest one up to now
		if (rate > biggestProductionRate)
			biggestProductionRate = rate;
	}

	// Loop on the combining list
	for (auto it = effCombiningList.begin(); it != effCombiningList.end();
			++it) {
		// Get the reactants
		combiningReactant = (*it).first;
		// Compute the reaction constant
		rate = calculateReactionRateConstant(*this, *combiningReactant);
		// Set it in the combining cluster
		(*it).kConstant = rate / (double) nTot;
	}

	// Loop on the dissociating list
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end();
			++it) {
		dissociatingCluster = (*it).first;
		// The second element of the pair is the cluster that is also
		// emitted by the dissociation
		otherEmittedCluster = (*it).second;
		// Compute the dissociation constant
		// The order of the cluster is important here because of the binding
		// energy used in the computation. It is taken from the type of the first cluster
		// which must be the single one
		// otherEmittedCluster is the single size one
		rate = calculateDissociationConstant(*dissociatingCluster,
				*otherEmittedCluster, *this);
		// Set it in the pair
		(*it).kConstant = rate / (double) nTot;
	}

	// Loop on the emission list
	for (auto it = effEmissionList.begin(); it != effEmissionList.end(); ++it) {
		firstCluster = (*it).first;
		secondCluster = (*it).second;
		// Compute the dissociation rate
		rate = calculateDissociationConstant(*this, *firstCluster,
				*secondCluster);
		// Set it in the pair
		(*it).kConstant = rate / (double) nTot;
	}

	// Set the biggest rate to the biggest production rate
	biggestRate = biggestProductionRate;

	return;
}

void PSISuperCluster::resetConnectivities() {
	// Clear both sets
	reactionConnectivitySet.clear();
	dissociationConnectivitySet.clear();

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);
	setDissociationConnectivity(id);
	setReactionConnectivity(heMomId);
	setDissociationConnectivity(heMomId);
	setReactionConnectivity(vMomId);
	setDissociationConnectivity(vMomId);

	// Loop over all the reacting pairs
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
		// The cluster is connecting to both clusters in the pair
		setReactionConnectivity((*it).first->getId());
		setReactionConnectivity((*it).first->getHeMomentumId());
		setReactionConnectivity((*it).first->getVMomentumId());
		setReactionConnectivity((*it).second->getId());
		setReactionConnectivity((*it).second->getHeMomentumId());
		setReactionConnectivity((*it).second->getVMomentumId());
	}

	// Loop over all the combining pairs
	for (auto it = effCombiningList.begin(); it != effCombiningList.end();
			++it) {
		// The cluster is connecting to the combining cluster
		setReactionConnectivity((*it).first->getId());
		setReactionConnectivity((*it).first->getHeMomentumId());
		setReactionConnectivity((*it).first->getVMomentumId());
	}

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end();
			++it) {
		// The cluster is connecting to the combining cluster
		setDissociationConnectivity((*it).first->getId());
		setDissociationConnectivity((*it).first->getHeMomentumId());
		setDissociationConnectivity((*it).first->getVMomentumId());
	}

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	// Initialize the partial vector for the momentum
	int dof = network->getDOF();
	heMomentumPartials.resize(dof, 0.0);
	vMomentumPartials.resize(dof, 0.0);

	return;
}

double PSISuperCluster::getTotalFlux() {
	// Initialize the fluxes
	heMomentumFlux = 0.0;
	vMomentumFlux = 0.0;

	// Get the fluxes
	double prodFlux = getProductionFlux();
	double dissFlux = getDissociationFlux();
	double combFlux = getCombinationFlux();
	double emissFlux = getEmissionFlux();

	return prodFlux - combFlux + dissFlux - emissFlux;
}

double PSISuperCluster::getDissociationFlux() {
	// Initial declarations
	double flux = 0.0, value = 0.0;
	PSICluster *dissociatingCluster = nullptr;

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end();
			++it) {
		// Get the dissociating clusters
		dissociatingCluster = (*it).first;
		double l0A = dissociatingCluster->getConcentration(0.0, 0.0);
		double lHeA = dissociatingCluster->getHeMomentum();
		double lVA = dissociatingCluster->getVMomentum();
		// Update the flux
		value = (*it).kConstant;
		flux += value * ((*it).a00 * l0A + (*it).a10 * lHeA + (*it).a20 * lVA);
		// Compute the momentum fluxes
		heMomentumFlux += value
				* ((*it).a01 * l0A + (*it).a11 * lHeA + (*it).a21 * lVA);
		vMomentumFlux += value
				* ((*it).a02 * l0A + (*it).a12 * lHeA + (*it).a22 * lVA);
	}

	// Return the flux
	return flux;
}

double PSISuperCluster::getEmissionFlux() {
	// Initial declarations
	double flux = 0.0, value = 0.0;

	// Loop over all the emission pairs
	for (auto it = effEmissionList.begin(); it != effEmissionList.end(); ++it) {
		// Update the flux
		value = (*it).kConstant;
		flux += value * ((*it).a00 * l0 + (*it).a10 * l1He + (*it).a20 * l1V);
		// Compute the momentum fluxes
		heMomentumFlux -= value
				* ((*it).a01 * l0 + (*it).a11 * l1He + (*it).a21 * l1V);
		vMomentumFlux -= value
				* ((*it).a02 * l0 + (*it).a12 * l1He + (*it).a22 * l1V);
	}

	return flux;
}

double PSISuperCluster::getProductionFlux() {
	// Local declarations
	double flux = 0.0, value = 0.0;
	PSICluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Loop over all the reacting pairs
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
		// Get the two reacting clusters
		firstReactant = (*it).first;
		secondReactant = (*it).second;
		double l0A = firstReactant->getConcentration(0.0, 0.0);
		double l0B = secondReactant->getConcentration(0.0, 0.0);
		double lHeA = firstReactant->getHeMomentum();
		double lHeB = secondReactant->getHeMomentum();
		double lVA = firstReactant->getVMomentum();
		double lVB = secondReactant->getVMomentum();
		// Update the flux
		value = (*it).kConstant;
		flux += value
				* ((*it).a000 * l0A * l0B + (*it).a010 * l0A * lHeB
						+ (*it).a020 * l0A * lVB + (*it).a100 * lHeA * l0B
						+ (*it).a110 * lHeA * lHeB + (*it).a120 * lHeA * lVB
						+ (*it).a200 * lVA * l0B + (*it).a210 * lVA * lHeB
						+ (*it).a220 * lVA * lVB);
		// Compute the momentum fluxes
		heMomentumFlux += value
				* ((*it).a001 * l0A * l0B + (*it).a011 * l0A * lHeB
						+ (*it).a021 * l0A * lVB + (*it).a101 * lHeA * l0B
						+ (*it).a111 * lHeA * lHeB + (*it).a121 * lHeA * lVB
						+ (*it).a201 * lVA * l0B + (*it).a211 * lVA * lHeB
						+ (*it).a221 * lVA * lVB);
		vMomentumFlux += value
				* ((*it).a002 * l0A * l0B + (*it).a012 * l0A * lHeB
						+ (*it).a022 * l0A * lVB + (*it).a102 * lHeA * l0B
						+ (*it).a112 * lHeA * lHeB + (*it).a122 * lHeA * lVB
						+ (*it).a202 * lVA * l0B + (*it).a212 * lVA * lHeB
						+ (*it).a222 * lVA * lVB);
	}

	// Return the production flux
	return flux;
}

double PSISuperCluster::getCombinationFlux() {
	// Local declarations
	double flux = 0.0, value = 0.0;
	PSICluster *combiningCluster = nullptr;

	// Loop over all the combining clusters
	for (auto it = effCombiningList.begin(); it != effCombiningList.end();
			++it) {
		// Get the two reacting clusters
		combiningCluster = (*it).first;
		double l0B = combiningCluster->getConcentration(0.0, 0.0);
		double lHeB = combiningCluster->getHeMomentum();
		double lVB = combiningCluster->getVMomentum();
		// Update the flux
		value = (*it).kConstant;
		flux += value
				* ((*it).a000 * l0B * l0 + (*it).a100 * l0B * l1He
						+ (*it).a200 * l0B * l1V + (*it).a010 * lHeB * l0
						+ (*it).a110 * lHeB * l1He + (*it).a210 * lHeB * l1V
						+ (*it).a020 * lVB * l0 + (*it).a120 * lVB * l1He
						+ (*it).a220 * lVB * l1V);
		// Compute the momentum fluxes
		heMomentumFlux -= value
				* ((*it).a001 * l0B * l0 + (*it).a101 * l0B * l1He
						+ (*it).a201 * l0B * l1V + (*it).a011 * lHeB * l0
						+ (*it).a111 * lHeB * l1He + (*it).a211 * lHeB * l1V
						+ (*it).a021 * lVB * l0 + (*it).a121 * lVB * l1He
						+ (*it).a221 * lVB * l1V);
		vMomentumFlux -= value
				* ((*it).a002 * l0B * l0 + (*it).a102 * l0B * l1He
						+ (*it).a202 * l0B * l1V + (*it).a012 * lHeB * l0
						+ (*it).a112 * lHeB * l1He + (*it).a212 * lHeB * l1V
						+ (*it).a022 * lVB * l0 + (*it).a122 * lVB * l1He
						+ (*it).a222 * lVB * l1V);
	}

	return flux;
}

void PSISuperCluster::getPartialDerivatives(
		std::vector<double> & partials) const {
	// Reinitialize the momentum partial derivatives vector
	std::fill(heMomentumPartials.begin(), heMomentumPartials.end(), 0.0);
	std::fill(vMomentumPartials.begin(), vMomentumPartials.end(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return;
}

void PSISuperCluster::getProductionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	double value = 0.0;
	int index = 0;
	PSICluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A

	// Loop over all the reacting pairs
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
		// Get the two reacting clusters
		firstReactant = (*it).first;
		secondReactant = (*it).second;
		double l0A = firstReactant->getConcentration(0.0, 0.0);
		double l0B = secondReactant->getConcentration(0.0, 0.0);
		double lHeA = firstReactant->getHeMomentum();
		double lHeB = secondReactant->getHeMomentum();
		double lVA = firstReactant->getVMomentum();
		double lVB = secondReactant->getVMomentum();

		// Compute the contribution from the first part of the reacting pair
		value = (*it).kConstant;
		index = firstReactant->getId() - 1;
		partials[index] += value
				* ((*it).a000 * l0B + (*it).a010 * lHeB + (*it).a020 * lVB);
		heMomentumPartials[index] += value
				* ((*it).a001 * l0B + (*it).a011 * lHeB + (*it).a021 * lVB);
		vMomentumPartials[index] += value
				* ((*it).a002 * l0B + (*it).a012 * lHeB + (*it).a022 * lVB);
		index = firstReactant->getHeMomentumId() - 1;
		partials[index] += value
				* ((*it).a100 * l0B + (*it).a110 * lHeB + (*it).a120 * lVB);
		heMomentumPartials[index] += value
				* ((*it).a101 * l0B + (*it).a111 * lHeB + (*it).a121 * lVB);
		vMomentumPartials[index] += value
				* ((*it).a102 * l0B + (*it).a112 * lHeB + (*it).a122 * lVB);
		index = firstReactant->getVMomentumId() - 1;
		partials[index] += value
				* ((*it).a200 * l0B + (*it).a210 * lHeB + (*it).a220 * lVB);
		heMomentumPartials[index] += value
				* ((*it).a201 * l0B + (*it).a211 * lHeB + (*it).a221 * lVB);
		vMomentumPartials[index] += value
				* ((*it).a202 * l0B + (*it).a212 * lHeB + (*it).a222 * lVB);
		// Compute the contribution from the second part of the reacting pair
		index = secondReactant->getId() - 1;
		partials[index] += value
				* ((*it).a000 * l0A + (*it).a100 * lHeA + (*it).a200 * lVA);
		heMomentumPartials[index] += value
				* ((*it).a001 * l0A + (*it).a101 * lHeA + (*it).a201 * lVA);
		vMomentumPartials[index] += value
				* ((*it).a002 * l0A + (*it).a102 * lHeA + (*it).a202 * lVA);
		index = secondReactant->getHeMomentumId() - 1;
		partials[index] += value
				* ((*it).a010 * l0A + (*it).a110 * lHeA + (*it).a210 * lVA);
		heMomentumPartials[index] += value
				* ((*it).a011 * l0A + (*it).a111 * lHeA + (*it).a211 * lVA);
		vMomentumPartials[index] += value
				* ((*it).a012 * l0A + (*it).a112 * lHeA + (*it).a212 * lVA);
		index = secondReactant->getVMomentumId() - 1;
		partials[index] += value
				* ((*it).a020 * l0A + (*it).a120 * lHeA + (*it).a220 * lVA);
		heMomentumPartials[index] += value
				* ((*it).a021 * l0A + (*it).a121 * lHeA + (*it).a221 * lVA);
		vMomentumPartials[index] += value
				* ((*it).a022 * l0A + (*it).a122 * lHeA + (*it).a222 * lVA);
	}

	return;
}

void PSISuperCluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int index = 0;
	PSICluster *cluster = nullptr;
	double value = 0.0;

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A

	// Loop over all the combining clusters
	for (auto it = effCombiningList.begin(); it != effCombiningList.end();
			++it) {
		// Get the two reacting clusters
		cluster = (*it).first;
		double l0B = cluster->getConcentration(0.0, 0.0);
		double lHeB = cluster->getHeMomentum();
		double lVB = cluster->getVMomentum();

		// Compute the contribution from the combining cluster
		value = (*it).kConstant;
		index = cluster->getId() - 1;
		partials[index] -= value
				* ((*it).a000 * l0 + (*it).a100 * l1He + (*it).a200 * l1V);
		heMomentumPartials[index] -= value
				* ((*it).a001 * l0 + (*it).a101 * l1He + (*it).a201 * l1V);
		vMomentumPartials[index] -= value
				* ((*it).a002 * l0 + (*it).a102 * l1He + (*it).a202 * l1V);
		index = cluster->getHeMomentumId() - 1;
		partials[index] -= value
				* ((*it).a010 * l0 + (*it).a110 * l1He + (*it).a210 * l1V);
		heMomentumPartials[index] -= value
				* ((*it).a011 * l0 + (*it).a111 * l1He + (*it).a211 * l1V);
		vMomentumPartials[index] -= value
				* ((*it).a012 * l0 + (*it).a112 * l1He + (*it).a212 * l1V);
		index = cluster->getVMomentumId() - 1;
		partials[index] -= value
				* ((*it).a020 * l0 + (*it).a120 * l1He + (*it).a220 * l1V);
		heMomentumPartials[index] -= value
				* ((*it).a021 * l0 + (*it).a121 * l1He + (*it).a221 * l1V);
		vMomentumPartials[index] -= value
				* ((*it).a022 * l0 + (*it).a122 * l1He + (*it).a222 * l1V);
		// Compute the contribution from this cluster
		index = id - 1;
		partials[index] -= value
				* ((*it).a000 * l0B + (*it).a010 * lHeB + (*it).a020 * lVB);
		heMomentumPartials[index] -= value
				* ((*it).a001 * l0B + (*it).a011 * lHeB + (*it).a021 * lVB);
		vMomentumPartials[index] -= value
				* ((*it).a002 * l0B + (*it).a012 * lHeB + (*it).a022 * lVB);
		index = heMomId - 1;
		partials[index] -= value
				* ((*it).a100 * l0B + (*it).a110 * lHeB + (*it).a120 * lVB);
		heMomentumPartials[index] -= value
				* ((*it).a101 * l0B + (*it).a111 * lHeB + (*it).a121 * lVB);
		vMomentumPartials[index] -= value
				* ((*it).a102 * l0B + (*it).a112 * lHeB + (*it).a122 * lVB);
		index = vMomId - 1;
		partials[index] -= value
				* ((*it).a200 * l0B + (*it).a210 * lHeB + (*it).a220 * lVB);
		heMomentumPartials[index] -= value
				* ((*it).a201 * l0B + (*it).a211 * lHeB + (*it).a221 * lVB);
		vMomentumPartials[index] -= value
				* ((*it).a202 * l0B + (*it).a212 * lHeB + (*it).a222 * lVB);
	}

	return;
}

void PSISuperCluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int index = 0;
	PSICluster *cluster = nullptr;
	double value = 0.0;

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end();
			++it) {
		// Get the dissociating clusters
		cluster = (*it).first;
		// Compute the contribution from the dissociating cluster
		value = (*it).kConstant;
		index = cluster->getId() - 1;
		partials[index] += value * ((*it).a00);
		heMomentumPartials[index] += value * ((*it).a01);
		vMomentumPartials[index] += value * ((*it).a02);
		index = cluster->getHeMomentumId() - 1;
		partials[index] += value * ((*it).a10);
		heMomentumPartials[index] += value * ((*it).a11);
		vMomentumPartials[index] += value * ((*it).a12);
		index = cluster->getVMomentumId() - 1;
		partials[index] += value * ((*it).a20);
		heMomentumPartials[index] += value * ((*it).a21);
		vMomentumPartials[index] += value * ((*it).a22);
	}

	return;
}

void PSISuperCluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int index = 0;
	double value = 0.0;

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)

	// Loop over all the emission pairs
	for (auto it = effEmissionList.begin(); it != effEmissionList.end(); ++it) {
		// Compute the contribution from the dissociating cluster
		value = (*it).kConstant;
		index = id - 1;
		partials[index] -= value * ((*it).a00);
		heMomentumPartials[index] -= value * ((*it).a01);
		vMomentumPartials[index] -= value * ((*it).a02);
		index = heMomId - 1;
		partials[index] -= value * ((*it).a10);
		heMomentumPartials[index] -= value * ((*it).a11);
		vMomentumPartials[index] -= value * ((*it).a12);
		index = vMomId - 1;
		partials[index] -= value * ((*it).a20);
		heMomentumPartials[index] -= value * ((*it).a21);
		vMomentumPartials[index] -= value * ((*it).a22);
	}

	return;
}

void PSISuperCluster::getHeMomentPartialDerivatives(
		std::vector<double> & partials) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = heMomentumPartials[i];
	}

	return;
}

void PSISuperCluster::getVMomentPartialDerivatives(
		std::vector<double> & partials) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = vMomentumPartials[i];
	}

	return;
}
