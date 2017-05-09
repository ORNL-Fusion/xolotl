// Includes
#include "AlloySuperCluster.h"
#include "AlloyClusterReactionNetwork.h"
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

/**
 * The moment partials.
 */
std::vector<double> momentPartials;

AlloySuperCluster::AlloySuperCluster(double num, int nTot, int width, double radius,
		double energy, std::string type, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		AlloyCluster(registry), numAtom(num), nTot(nTot), l0(0.0), l1(0.0), dispersion(
				0.0), momentFlux(0.0) {
	// Set the cluster size
	size = (int) numAtom;

	// Update the composition map
	compositionMap[type] = (int) (numAtom * (double) nTot);

	// Set the width
	sectionWidth = width;

	// Set the reaction radius and formation energy
	reactionRadius = radius;
	formationEnergy = energy;

	// Set the diffusion factor and the migration energy
	migrationEnergy = std::numeric_limits<double>::infinity();
	diffusionFactor = 0.0;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << type << "_" << numAtom;
	name = nameStream.str();
	// Set the typename appropriately
	if (type == voidType) typeName = AlloyVoidSuperType;
	else if (type == faultedType) typeName = AlloyFaultedSuperType;
	else if (type == frankType) typeName = AlloyFrankSuperType;
	else if (type == perfectType) typeName = AlloyPerfectSuperType;

	return;
}

AlloySuperCluster::AlloySuperCluster(AlloySuperCluster &other) :
		AlloyCluster(other) {
	numAtom = other.numAtom;
	nTot = other.nTot;
	sectionWidth = other.sectionWidth;
	l0 = other.l0;
	l1 = other.l1;
	dispersion = other.dispersion;
	reactingMap = other.reactingMap;
	combiningMap = other.combiningMap;
	dissociatingMap = other.dissociatingMap;
	emissionMap = other.emissionMap;
	effReactingList = other.effReactingList;
	effCombiningList = other.effCombiningList;
	effDissociatingList = other.effDissociatingList;
	effEmissionList = other.effEmissionList;
	momentFlux = other.momentFlux;

	return;
}

std::shared_ptr<IReactant> AlloySuperCluster::clone() {
	std::shared_ptr<IReactant> reactant(new AlloySuperCluster(*this));

	return reactant;
}

void AlloySuperCluster::setReactionNetwork(
		const std::shared_ptr<IReactionNetwork> reactionNetwork) {
	// Call the superclass's method to actually set the reference
	Reactant::setReactionNetwork(reactionNetwork);

	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	// Aggregate the reacting pairs and combining reactants from the atomVector
	// Loop on the atomVector
	for (int i = 0; i < atomVector.size(); i++) {
		// Get the cluster size
		double size = atomVector[i]->getSize();
		// Get all vectors
		auto react = atomVector[i]->reactingPairs;
		auto combi = atomVector[i]->combiningReactants;
		auto disso = atomVector[i]->dissociatingPairs;
		auto emi = atomVector[i]->emissionPairs;

		// Set them in the super cluster map
		reactingMap[size] = react;
		combiningMap[size] = combi;
		dissociatingMap[size] = disso;
		emissionMap[size] = emi;
	}

	// Compute the dispersion
	computeDispersion();

	return;
}

double AlloySuperCluster::getConcentration(double distAtom, double distB) const {
	return l0 + (distAtom * l1);
}

double AlloySuperCluster::getMoment() const {
	return l1;
}

double AlloySuperCluster::getTotalConcentration() const {
	// Initial declarations
	int index = 0;
	double distance = 0.0, conc = 0.0;

	// Loop on the xenon width
	for (int k = 0; k < sectionWidth; k++) {
		// Compute the xenon index
		index = (int) (numAtom - (double) sectionWidth / 2.0) + k + 1;

		// Check if this cluster exists
		if (reactingMap.find(index) == reactingMap.end())
			continue;

		// Compute the distances
		distance = getDistance(index);

		// Add the concentration of each cluster in the group
		conc += getConcentration(distance);
	}

	return conc;
}

double AlloySuperCluster::getTotalAtomConcentration() const {
	// Initial declarations
	int index = 0;
	double distance = 0.0, conc = 0.0;

	// Loop on the xenon width
	for (int k = 0; k < sectionWidth; k++) {
		// Compute the xenon index
		index = (int) (numAtom - (double) sectionWidth / 2.0) + k + 1;

		// Check if this cluster exists
		if (reactingMap.find(index) == reactingMap.end())
			continue;

		// Compute the distances
		distance = getDistance(index);

		// Add the concentration of each cluster in the group times the number of xenon atoms
		conc += getConcentration(distance) * (double) index;
	}

	return conc;
}

double AlloySuperCluster::getDistance(int xe) const {
	if (sectionWidth == 1)
		return 0.0;
	return 2.0 * (double) (xe - numAtom) / ((double) sectionWidth - 1.0);
}

void AlloySuperCluster::computeDispersion() {
	// Local declarations
	AlloyCluster *firstReactant, *secondReactant, *combiningReactant,
			*dissociatingCluster, *otherEmittedCluster, *firstCluster,
			*secondCluster;
	double rate = 0.0;
	int index = 0;
	// Initialize the dispersion sum
	double nAtomSquare = 0.0;

	// Loop on the xenon width
	for (int k = 0; k < sectionWidth; k++) {
		// Compute the vacancy index
		index = (int) (numAtom - (double) sectionWidth / 2.0) + k + 1;

		// Check if this cluster exists
		if (reactingMap.find(index) == reactingMap.end())
			continue;

		// Compute nSquare for the dispersion
		nAtomSquare += (double) index * index;
	}

	// Compute the dispersion
	if (sectionWidth == 1)
		dispersion = 1.0;
	else {
		dispersion = 2.0
				* (nAtomSquare
						- ((double) size * (double) nTot
								* ((double) size * (double) nTot
										/ (double) sectionWidth)))
				/ ((double) (sectionWidth * (sectionWidth - 1)));
	}

	return;
}

void AlloySuperCluster::optimizeReactions() {
	// Local declarations
	double factor = 0.0, distance = 0.0;
	AlloyCluster *firstReactant, *secondReactant, *combiningReactant,
			*dissociatingCluster, *otherEmittedCluster, *firstCluster,
			*secondCluster;
	int index = 0;

	// Loop on the reacting map
	for (auto mapIt = reactingMap.begin(); mapIt != reactingMap.end();
			++mapIt) {
		// Get the pairs
		auto pairs = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = pairs.begin(); it != pairs.end();) {
			// Get the two reacting clusters
			firstReactant = (*it).first;
			secondReactant = (*it).second;

			// Create the corresponding production reaction
			auto reaction = std::make_shared<ProductionReaction>(firstReactant,
					secondReactant);
			// Add it to the network
			reaction = network->addProductionReaction(reaction);

			// Create a new SuperClusterProductionPair
			SuperClusterProductionPair superPair(firstReactant, secondReactant,
					reaction.get());

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != reactingMap.end();
					++mapItBis) {
				// Compute the xenon index
				index = mapItBis->first;
				factor = (double) (index - numAtom) / dispersion;

				// Get the pairs
				auto pairsBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = pairsBis.begin(); itBis != pairsBis.end();) {
					// Get the two reacting clusters
					auto firstReactantBis = (*itBis).first;
					auto secondReactantBis = (*itBis).second;

					// Check if it is the same reaction
					if (firstReactantBis == firstReactant
							&& secondReactantBis == secondReactant) {
						superPair.a000 += 1.0;
						superPair.a001 += factor;
						superPair.a100 += (*itBis).firstDistance;
						superPair.a101 += (*itBis).firstDistance * factor;
						superPair.a010 += (*itBis).secondDistance;
						superPair.a011 += (*itBis).secondDistance * factor;
						superPair.a110 += (*itBis).firstDistance
								* (*itBis).secondDistance;
						superPair.a111 += (*itBis).firstDistance
								* (*itBis).secondDistance * factor;

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

	// Loop on the combining map
	for (auto mapIt = combiningMap.begin(); mapIt != combiningMap.end();
			++mapIt) {
		// Get the pairs
		auto clusters = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = clusters.begin(); it != clusters.end();) {
			// Get the combining cluster
			combiningReactant = (*it).combining;

			// Create the corresponding production reaction
			auto reaction = std::make_shared<ProductionReaction>(this,
					combiningReactant);
			// Add it to the network
			reaction = network->addProductionReaction(reaction);

			// Create a new SuperClusterProductionPair with NULL as the second cluster because
			// we do not need it
			SuperClusterProductionPair superPair(combiningReactant, NULL,
					reaction.get());

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != combiningMap.end();
					++mapItBis) {
				// Compute the xenon index
				index = mapItBis->first;
				distance = getDistance(index);
				factor = (double) (index - numAtom) / dispersion;

				// Get the pairs
				auto clustersBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = clustersBis.begin();
						itBis != clustersBis.end();) {
					// Get the two reacting clusters
					auto combiningReactantBis = (*itBis).combining;

					// Check if it is the same reaction
					if (combiningReactantBis == combiningReactant) {
						superPair.a000 += 1.0;
						superPair.a001 += factor;
						superPair.a010 += (*itBis).distance;
						superPair.a011 += (*itBis).distance * factor;
						superPair.a100 += distance;
						superPair.a101 += distance * factor;
						superPair.a110 += (*itBis).distance * distance;
						superPair.a111 += (*itBis).distance * distance * factor;

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

	// Loop on the dissociating map
	for (auto mapIt = dissociatingMap.begin(); mapIt != dissociatingMap.end();
			++mapIt) {
		// Get the pairs
		auto pairs = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = pairs.begin(); it != pairs.end();) {
			// Get the two reacting clusters
			dissociatingCluster = (*it).first;
			otherEmittedCluster = (*it).second;

			// Create a dissociation reaction
			auto reaction = std::make_shared<DissociationReaction>(
					dissociatingCluster, this, otherEmittedCluster);
			// Add it to the network
			reaction = network->addDissociationReaction(reaction);

			// Create a new SuperClusterProductionPair
			SuperClusterDissociationPair superPair(dissociatingCluster,
					otherEmittedCluster, reaction.get());

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != dissociatingMap.end();
					++mapItBis) {
				// Compute the xenon index
				index = mapItBis->first;
				factor = (double) (index - numAtom) / dispersion;

				// Get the pairs
				auto pairsBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = pairsBis.begin(); itBis != pairsBis.end();) {
					// Get the two reacting clusters
					auto dissociatingClusterBis = (*itBis).first;
					auto otherEmittedClusterBis = (*itBis).second;

					// Check if it is the same reaction
					if (dissociatingClusterBis == dissociatingCluster
							&& otherEmittedClusterBis == otherEmittedCluster) {
						superPair.a00 += 1.0;
						superPair.a01 += factor;
						superPair.a10 += (*itBis).firstDistance;
						superPair.a11 += (*itBis).firstDistance * factor;

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

	// Loop on the emission map
	for (auto mapIt = emissionMap.begin(); mapIt != emissionMap.end();
			++mapIt) {
		// Get the pairs
		auto pairs = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = pairs.begin(); it != pairs.end();) {
			// Get the two reacting clusters
			firstCluster = (*it).first;
			secondCluster = (*it).second;

			// Create a dissociation reaction
			auto reaction = std::make_shared<DissociationReaction>(this,
					firstCluster, secondCluster);
			// Add it to the network
			reaction = network->addDissociationReaction(reaction);

			// Create a new SuperClusterProductionPair
			SuperClusterDissociationPair superPair(firstCluster, secondCluster,
					reaction.get());

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != emissionMap.end();
					++mapItBis) {
				// Compute the xenon index
				index = mapItBis->first;
				distance = getDistance(index);
				factor = (double) (index - numAtom) / dispersion;

				// Get the pairs
				auto pairsBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = pairsBis.begin(); itBis != pairsBis.end();) {
					// Get the two reacting clusters
					auto firstClusterBis = (*itBis).first;
					auto secondClusterBis = (*itBis).second;

					// Check if it is the same reaction
					if (firstClusterBis == firstCluster
							&& secondClusterBis == secondCluster) {
						superPair.a00 += 1.0;
						superPair.a01 += factor;
						superPair.a10 += distance;
						superPair.a11 += distance * factor;

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
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	return;
}

void AlloySuperCluster::resetConnectivities() {
	// Clear both sets
	reactionConnectivitySet.clear();
	dissociationConnectivitySet.clear();

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);
	setDissociationConnectivity(id);
	setReactionConnectivity(momId);
	setDissociationConnectivity(momId);

	// Loop over all the reacting pairs
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
		// The cluster is connecting to both clusters in the pair
		setReactionConnectivity((*it).first->getId());
		setReactionConnectivity((*it).first->getMomentId());
		setReactionConnectivity((*it).second->getId());
		setReactionConnectivity((*it).second->getMomentId());
	}

	// Loop over all the combining pairs
	for (auto it = effCombiningList.begin(); it != effCombiningList.end();
			++it) {
		// The cluster is connecting to the combining cluster
		setReactionConnectivity((*it).first->getId());
		setReactionConnectivity((*it).first->getMomentId());
	}

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end();
			++it) {
		// The cluster is connecting to the combining cluster
		setDissociationConnectivity((*it).first->getId());
		setDissociationConnectivity((*it).first->getMomentId());
	}

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	// Initialize the partial vector for the moment
	int dof = network->getDOF();
	momentPartials.resize(dof, 0.0);

	return;
}

double AlloySuperCluster::getTotalFlux() {
	// Initialize the moment flux
	momentFlux = 0.0;

	// Get the fluxes
	double prodFlux = getProductionFlux();
	double dissFlux = getDissociationFlux();
	double combFlux = getCombinationFlux();
	double emissFlux = getEmissionFlux();

	return prodFlux - combFlux + dissFlux - emissFlux;
}

double AlloySuperCluster::getDissociationFlux() {
	// Initial declarations
	double flux = 0.0, value = 0.0;
	AlloyCluster *dissociatingCluster = nullptr;

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end();
			++it) {
		// Get the dissociating clusters
		dissociatingCluster = (*it).first;
		double l0A = dissociatingCluster->getConcentration(0.0);
		double l1A = dissociatingCluster->getMoment();
		// Update the flux
		value = *((*it).kConstant) / (double) nTot;
		flux += value * ((*it).a00 * l0A + (*it).a10 * l1A);
		// Compute the moment fluxes
		momentFlux += value * ((*it).a01 * l0A + (*it).a11 * l1A);
	}

	// Return the flux
	return flux;
}

double AlloySuperCluster::getEmissionFlux() {
	// Initial declarations
	double flux = 0.0, value = 0.0;

	// Loop over all the emission pairs
	for (auto it = effEmissionList.begin(); it != effEmissionList.end(); ++it) {
		// Update the flux
		value = *((*it).kConstant) / (double) nTot;
		flux += value * ((*it).a00 * l0 + (*it).a10 * l1);
		// Compute the moment fluxes
		momentFlux -= value * ((*it).a01 * l0 + (*it).a11 * l1);
	}

	return flux;
}

double AlloySuperCluster::getProductionFlux() {
	// Local declarations
	double flux = 0.0, value = 0.0;
	AlloyCluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Loop over all the reacting pairs
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
		// Get the two reacting clusters
		firstReactant = (*it).first;
		secondReactant = (*it).second;
		double l0A = firstReactant->getConcentration();
		double l0B = secondReactant->getConcentration();
		double l1A = firstReactant->getMoment();
		double l1B = secondReactant->getMoment();
		// Update the flux
		value = *((*it).kConstant) / (double) nTot;
		flux += value
				* ((*it).a000 * l0A * l0B + (*it).a010 * l0A * l1B
						+ (*it).a100 * l1A * l0B + (*it).a110 * l1A);
		// Compute the moment flux
		momentFlux += value
				* ((*it).a001 * l0A * l0B + (*it).a011 * l0A * l1B
						+ (*it).a101 * l1A * l0B + (*it).a111 * l1A);
	}

	// Return the production flux
	return flux;
}

double AlloySuperCluster::getCombinationFlux() {
	// Local declarations
	double flux = 0.0, value = 0.0;
	AlloyCluster *combiningCluster = nullptr;

	// Loop over all the combining clusters
	for (auto it = effCombiningList.begin(); it != effCombiningList.end();
			++it) {
		// Get the two reacting clusters
		combiningCluster = (*it).first;
		double l0A = combiningCluster->getConcentration();
		double l1A = combiningCluster->getMoment();
		// Update the flux
		value = *((*it).kConstant) / (double) nTot;
		flux += value
				* ((*it).a000 * l0A * l0 + (*it).a100 * l0A * l1
						+ (*it).a010 * l1A * l0 + (*it).a110 * l1A * l1);
		// Compute the moment flux
		momentFlux -= value
				* ((*it).a001 * l0A * l0 + (*it).a101 * l0A * l1
						+ (*it).a011 * l1A * l0 + (*it).a111 * l1A * l1);
	}

	return flux;
}

void AlloySuperCluster::getPartialDerivatives(
		std::vector<double> & partials) const {
	// Reinitialize the moment partial derivatives vector
	std::fill(momentPartials.begin(), momentPartials.end(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return;
}

void AlloySuperCluster::getProductionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	double value = 0.0;
	int index = 0;
	AlloyCluster *firstReactant = nullptr, *secondReactant = nullptr;

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
		double l0A = firstReactant->getConcentration();
		double l0B = secondReactant->getConcentration();
		double l1A = firstReactant->getMoment();
		double l1B = secondReactant->getMoment();

		// Compute the contribution from the first part of the reacting pair
		value = *((*it).kConstant) / (double) nTot;

//		std::cout << name << " : " << firstReactant->getName() << " + " << secondReactant->getName() << " " << l0A << " " << l0B << std::endl;

		index = firstReactant->getId() - 1;
		partials[index] += value * ((*it).a000 * l0B + (*it).a010 * l1B);
		momentPartials[index] += value
				* ((*it).a001 * l0B + (*it).a011 * l1B);
		index = firstReactant->getMomentId() - 1;
		partials[index] += value * ((*it).a100 * l0B + (*it).a110 * l1B);
		momentPartials[index] += value
				* ((*it).a101 * l0B + (*it).a111 * l1B);
		// Compute the contribution from the second part of the reacting pair
		index = secondReactant->getId() - 1;
		partials[index] += value * ((*it).a000 * l0A + (*it).a100 * l1A);
		momentPartials[index] += value
				* ((*it).a001 * l0A + (*it).a101 * l1A);
		index = secondReactant->getMomentId() - 1;
		partials[index] += value * ((*it).a010 * l0A + (*it).a110 * l1A);
		momentPartials[index] += value
				* ((*it).a011 * l0A + (*it).a111 * l1A);
	}

	return;
}

void AlloySuperCluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int index = 0;
	AlloyCluster *cluster = nullptr;
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
		double l0A = cluster->getConcentration();
		double l1A = cluster->getMoment();

		// Compute the contribution from the combining cluster
		value = *((*it).kConstant) / (double) nTot;
		index = cluster->getId() - 1;
		partials[index] -= value * ((*it).a000 * l0 + (*it).a100 * l1);
		momentPartials[index] -= value * ((*it).a001 * l0 + (*it).a101 * l1);
		index = cluster->getMomentId() - 1;
		partials[index] -= value * ((*it).a010 * l0 + (*it).a110 * l1);
		momentPartials[index] -= value * ((*it).a011 * l0 + (*it).a111 * l1);
		// Compute the contribution from this cluster
		index = id - 1;
		partials[index] -= value * ((*it).a000 * l0A + (*it).a010 * l1A);
		momentPartials[index] -= value
				* ((*it).a001 * l0A + (*it).a011 * l1A);
		index = momId - 1;
		partials[index] -= value * ((*it).a100 * l0A + (*it).a110 * l1A);
		momentPartials[index] -= value
				* ((*it).a101 * l0A + (*it).a111 * l1A);
	}

	return;
}

void AlloySuperCluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int index = 0;
	AlloyCluster *cluster = nullptr;
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
		value = *((*it).kConstant) / (double) nTot;
		index = cluster->getId() - 1;
		partials[index] += value * ((*it).a00);
		momentPartials[index] += value * ((*it).a01);
		index = cluster->getMomentId() - 1;
		partials[index] += value * ((*it).a10);
		momentPartials[index] += value * ((*it).a11);
	}

	return;
}

void AlloySuperCluster::getEmissionPartialDerivatives(
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
		value = *((*it).kConstant) / (double) nTot;
		index = id - 1;
		partials[index] -= value * ((*it).a00);
		momentPartials[index] -= value * ((*it).a01);
		index = momId - 1;
		partials[index] -= value * ((*it).a10);
		momentPartials[index] -= value * ((*it).a11);
	}

	return;
}

void AlloySuperCluster::getMomentPartialDerivatives(
		std::vector<double> & partials) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = momentPartials[i];
	}

	return;
}
