#include "PSICluster.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

PSICluster::PSICluster() {
	// Set the size
	size = 1;
	// Zero out the binding energies
	bindingEnergies.resize(4);
	bindingEnergies[0] = 0.0;
	bindingEnergies[1] = 0.0;
	bindingEnergies[2] = 0.0;
	bindingEnergies[3] = 0.0;
	// Zero out the diffusion factor and migration energy
	diffusionFactor = 0.0;
	migrationEnergy = 0.0;
	// Set the reactant name appropriately
	name = "PSICluster";
}

PSICluster::PSICluster(const int clusterSize) :
		Reactant() {

	// Set the size
	size = (clusterSize > 0) ? clusterSize : 1;
	// Zero out the binding energies
	bindingEnergies.resize(4);
	bindingEnergies[0] = 0.0;
	bindingEnergies[1] = 0.0;
	bindingEnergies[2] = 0.0;
	bindingEnergies[3] = 0.0;
	// Zero out the diffusion factor and migration energy
	diffusionFactor = 0.0;
	migrationEnergy = 0.0;
	// Set the reactant name appropriately
	name = "PSICluster";
}

// The copy constructor with a huge initialization list!
PSICluster::PSICluster(const PSICluster &other) :
		Reactant(other), size(other.size), diffusionFactor(
				other.diffusionFactor), bindingEnergies(other.bindingEnergies), migrationEnergy(
				other.migrationEnergy), reactionConnectivity(
				other.reactionConnectivity), dissociationConnectivity(
				other.dissociationConnectivity), reactingPairs(
				other.reactingPairs), combiningReactants(
				other.combiningReactants) {
}

std::shared_ptr<Reactant> PSICluster::clone() {
	std::shared_ptr<Reactant> reactant(new PSICluster(*this));
	return reactant;
}

PSICluster::~PSICluster() {
}

int PSICluster::getSize() const {
	// Return this cluster's size
	return size;
}

void PSICluster::setReactionNetwork(
		const std::shared_ptr<ReactionNetwork> reactionNetwork) {

	// Call the superclass's method to actually set the reference
	Reactant::setReactionNetwork(reactionNetwork);

	// ----- Handle the connectivty for PSIClusters -----

	// Extract properties from the network
	auto properties = network->getProperties();
	int connectivityLength = network->size();

	// Get the enabled reaction type flags
	bool reactionsEnabled = (properties["reactionsEnabled"] == "true");
	bool dissociationsEnabled = (properties["dissociationsEnabled"] == "true");

	// Resize the connectivity arrays
	reactionConnectivity.resize(connectivityLength, 0);
	dissociationConnectivity.resize(connectivityLength, 0);
	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();

	// Generate the reactant and dissociation connectivity arrays.
	// This only must be done once since the arrays are stored as
	// member attributes. Only perform these tasks if the reaction
	// types are enabled.
	if (reactionsEnabled) {
		createReactionConnectivity();
	}
	if (dissociationsEnabled) {
		createDissociationConnectivity();
	}
}

double PSICluster::getDissociationFlux(double temperature) const {

	int nReactants = 0, oneIndex = -1;
	double diss = 0.0, conc = 0.0;
	std::shared_ptr<PSICluster> first, second;
	std::map<std::string, int> oneHe, oneV, oneI;

	// Only try this if the network is available
	if (network != NULL) {
		// Get the total number of reactants in the network
		nReactants = network->size();
		// Get this PSICluster or subclasses' cluster map
		auto composition = getComposition();
		// Get the number of species to determine if this
		// cluster is mixed or single
		int numSpecies = (composition["He"] > 0) + (composition["V"] > 0)
				+ (composition["I"] > 0);
		// If no species, throw error
		if (numSpecies == 0) {
			// Bad if we have no species
			throw std::string("Cluster map contains no species");
		} else if (numSpecies == 1) {
			// We know we are a single species,
			// but we need to know which one so we can
			// get the correct species He, V, or I to calculate
			// the dissociation constant.
			if (composition["He"]) {
				second = std::dynamic_pointer_cast<PSICluster>(
						network->get("He", 1));
			} else if (composition["V"]) {
				second = std::dynamic_pointer_cast<PSICluster>(
						network->get("V", 1));
			} else if (composition["I"]) {
				second = std::dynamic_pointer_cast<PSICluster>(
						network->get("I", 1));
			}
			// Loop over all reactants and see if we
			// have a dissociation connection
			auto reactants = network->getAll(); //Assumed it is ordered exactly!
			for (int i = 0; i < nReactants; i++) {
				// Only calculate if we are connected
				if (dissociationConnectivity.at(i) == 1) {
					first = std::dynamic_pointer_cast<PSICluster>(
							reactants->at(i));
					// Calculate the dissociation flux
					diss += calculateDissociationConstant(*first, *second,
							temperature) * first->getConcentration();
				}
			}
		} else if (numSpecies == 2) {
			std::cout << "PSICluster Message: "
					<< "Caught invalid single-species composition! "
					<< composition["He"] << " " << composition["V"] << " "
					<< composition["I"] << std::endl;
			throw std::string(
					"Mixed Species dissociation flux must be implemented by subclass.");
		}

	}

	// Return the flux
	return diss;
}

double PSICluster::getProductionFlux(double temperature) const {

	// Local declarations
	double flux = 0.0, kPlus = 0.0;
	int thisClusterIndex = 0;
	std::shared_ptr<PSICluster> firstReactant, secondReactant;
	int nPairs = 0;

	// Only try this if the network is available
	if (network != NULL) {
		// Set the total number of reacting pairs
		nPairs = reactingPairs.size();
		// Loop over all the reacting pairs
		for (int i = 0; i < nPairs; i++) {
			// Get the reactants
			firstReactant = reactingPairs.at(i).first;
			secondReactant = reactingPairs.at(i).second;
			// Update the flux
			flux += calculateReactionRateConstant(*firstReactant,
					*secondReactant, temperature)
					* firstReactant->getConcentration()
					* secondReactant->getConcentration();
		}

	}

	// Return the production flux
	return flux;
}

/**
 * This operation returns the total change in this cluster due to
 * the combination of this cluster with others.
 * @param temperature The temperature at which to calculate the flux
 * @return The flux due to this cluster combining with other clusters.
 */
double PSICluster::getCombinationFlux(double temperature) const {
	// Local declarations
	double flux = 0.0, kPlus = 0.0, conc = 0.0;
	int thisClusterIndex = 0;
	std::shared_ptr<PSICluster> outerReactant, thisReactant;
	int nReactants = 0;

	// Only try this if the network is available
	if (network != NULL) {
		// Set the total network nReactants
		nReactants = combiningReactants.size();
		// Loop over all possible clusters
		auto composition = getComposition();
//		std::cout << "PSICluster Message 1: " << composition["He"] << " "
//				<< composition["V"] << " " << composition["I"] << std::endl;
		for (int j = 0; j < nReactants; j++) {
			outerReactant = std::dynamic_pointer_cast<PSICluster>(
					combiningReactants.at(j));
			conc = outerReactant->getConcentration();
			composition = outerReactant->getComposition();
//			std::cout << "PSICluster Message 2: " << composition["He"] << " "
//					<< composition["V"] << " " << composition["I"] << std::endl;
			// Calculate Second term of production flux
			flux += calculateReactionRateConstant(*this, *outerReactant,
					temperature) * conc;
		}

	}

	// Return the production flux
	return (flux * getConcentration());
}

double PSICluster::getTotalFlux(double temperature) const {
	return this->getProductionFlux(temperature)
			- this->getCombinationFlux(temperature)
			+ this->getDissociationFlux(temperature);
}

double PSICluster::getDiffusionFactor() const {
	// Return the diffusion factor
	return diffusionFactor;
}

void PSICluster::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;
	return;
}

double PSICluster::getDiffusionCoefficient(double temperature) const {
	// Use the Arrhenius equation to compute the diffusion coefficient
	double k_b = xolotlCore::kBoltzmann;
	double kernel = -migrationEnergy / (k_b * temperature);
	return diffusionFactor * exp(kernel);
}

std::vector<double> PSICluster::getBindingEnergies() const {
	// Local Declarations
	std::vector<double> energyVector;

	// Set the return binding energies
	energyVector = bindingEnergies;

	// Return the energies
	return energyVector;
}

void PSICluster::setBindingEnergies(const std::vector<double> energies) {
	// Set the binding energies
	bindingEnergies = energies;
	return;
}

double PSICluster::getMigrationEnergy() const {
	// Return the migration energy
	return migrationEnergy;
}

void PSICluster::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;
	return;
}

double PSICluster::calculateReactionRateConstant(
		const PSICluster & firstReactant, const PSICluster & secondReactant,
		const double temperature) const {

	// Get the reaction radii
	double r_first = firstReactant.getReactionRadius();
	double r_second = secondReactant.getReactionRadius();

	// Get the diffusion coefficients
	double firstDiffusion = firstReactant.getDiffusionCoefficient(temperature);
	double secondDiffusion = secondReactant.getDiffusionCoefficient(
			temperature);

	// Calculate and return
	double k_plus = 4.0 * xolotlCore::pi * (r_first + r_second)
			* (firstDiffusion + secondDiffusion);
	return k_plus;
}

double PSICluster::calculateDissociationConstant(
		const PSICluster & firstReactant, const PSICluster & secondReactant,
		const double temperature) const {

	// Local Declarations
	int bindingEnergyIndex = -1;
	double ra = 1, rb = 1;
	double atomicVolume = 1.0;
	auto secondComposition = secondReactant.getComposition();

	// Get the binding energy index
	if (secondComposition["He"] == 1 && secondComposition["V"] == 0
			&& secondComposition["I"] == 0) {
		bindingEnergyIndex = 0;
	} else if (secondComposition["He"] == 0 && secondComposition["V"] == 1
			&& secondComposition["I"] == 0) {
		bindingEnergyIndex = 1;
	} else if (secondComposition["He"] == 0 && secondComposition["V"] == 0
			&& secondComposition["I"] == 1) {
		bindingEnergyIndex = 2;
	} else {
		return 0.0;
	}

	// Calculate the Reaction Rate Constant
	double kPlus = calculateReactionRateConstant(firstReactant, secondReactant,
			temperature);

	// Calculate and return
	double k_minus_exp = exp(
			-1.0 * bindingEnergies.at(bindingEnergyIndex)
					/ (xolotlCore::kBoltzmann * temperature));
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;
	return k_minus;
}

bool PSICluster::isProductReactant(const Reactant & reactantI,
		const Reactant & reactantJ) {
	// Base class should just return false
	return false;
}

double PSICluster::getReactionRadius() const {
	return 0.001; // FIXME!
}

std::vector<int> PSICluster::getConnectivity() const {

	int connectivityLength = network->size();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);

	// The reaction and dissociation vectors must have a length equal to the
	// number of clusters
	if (reactionConnectivity.size() != connectivityLength) {
		throw std::string("The reaction vector is an incorrect length");
	}
	if (dissociationConnectivity.size() != connectivityLength) {
		throw std::string("The dissociation vector is an incorrect length");
	}

	// Merge the two vectors such that the final vector contains
	// a 1 at a position if either of the connectivity arrays
	// have a 1
	for (int i = 0; i < connectivityLength; i++) {
		// Consider each connectivity array only if its type is enabled
		connectivity[i] = reactionConnectivity[i]
				|| dissociationConnectivity[i];
	}

	// DEBUG
	std::cout << std::endl << "Reactant Column: " << this->name << "_"
			<< this->size << std::endl;

	// DEBUG
	for (int i = 0; i < connectivity.size(); i++) {
		std::cout << connectivity[i] << " ";
	}
	std::cout << std::endl;

	return connectivity;
}

void PSICluster::createReactionConnectivity() {
	// By default, generate an array with a zero for each reactant
	// in the network
	reactionConnectivity.clear();
	reactionConnectivity.resize(network->size(), 0);
}

void PSICluster::createDissociationConnectivity() {

	// Local Declarations
	int index = 0;
	std::shared_ptr<Reactant> smallerReactant, singleReactant;

	// ----- X_a --> X_(a-1) + X ------
	smallerReactant = network->get(name, size - 1);
	singleReactant = network->get(name, 1);
	dissociateClusters(singleReactant,smallerReactant);

	return;
}

/**
 * This operation creates the two dissociated clusters from this cluster.
 * It is called by createDissociationConnectivity to process the reaction
 * and handle the connectivity.
 * @param firstDissociatedCluster The first cluster removed by
 * dissociation.
 * @param secondDissociatedCluster The second cluster removed by
 * dissociation.
 */
void PSICluster::dissociateClusters(const std::shared_ptr<Reactant> & firstDissociatedCluster,
		const std::shared_ptr<Reactant> & secondDissociatedCluster) {

	int index = 0;

	if (firstDissociatedCluster && secondDissociatedCluster) {
		// Add the two reactants to the set. "He" has very simple dissociation
		// rules.
		index = network->getReactantId(*firstDissociatedCluster) - 1;
		dissociationConnectivity[index] = 1;
		index = network->getReactantId(*secondDissociatedCluster) - 1;
		dissociationConnectivity[index] = 1;
		// Connect this cluster to itself since any reaction will affect it
		index = network->getReactantId(*this) - 1;
		dissociationConnectivity[index] = 1;
	}

	return;
}

/**
 * This operation returns the list of partial derivatives of this Reactant
 * with respect to all other reactants in the network. The combined lists
 * of partial derivatives from all of the reactants in the network can be
 * used to form, for example, a Jacobian.
 *
 * @param the temperature at which the reactions are occurring
 * @return The partial derivatives for this reactant where index zero
 * corresponds to the first reactant in the list returned by the
 * ReactionNetwork::getAll() operation.
 */
std::vector<double> PSICluster::getPartialDerivatives(
		double temperature) const {

	// Create the array and fill it with zeros
	int length = network->size(), size = 0, index = 0;
	std::vector<double> partialDerivatives = std::vector<double>(length, 0.0);
	ReactingPair pair;
	std::shared_ptr<PSICluster> cluster;

	// Load up everything from the reacting pairs array
	size = reactingPairs.size();
	for (int i = 0; i < size; i++) {
		pair = reactingPairs[i];
		// Compute the contribution from the first part of the reacting pair
		index = network->getReactantId(*(pair.first)) - 1;
		partialDerivatives[index] += calculateReactionRateConstant(*this,
				*(pair.first), temperature);
		// Compute the contribution from the second part of the reacting pair
		index = network->getReactantId(*(pair.second)) - 1;
		partialDerivatives[index] += calculateReactionRateConstant(*this,
				*(pair.second), temperature);
		std::cout << "Partial Derivative = " << partialDerivatives[index] << std::endl;
	}

	// Load up everything from the combining reactants
	size = combiningReactants.size();
	for (int i = 0; i < size; i++) {
		cluster = std::dynamic_pointer_cast<PSICluster>(combiningReactants[i]);
		// Compute the contribution from the cluster
		index = network->getReactantId(*cluster) - 1;
		partialDerivatives[index] += calculateReactionRateConstant(*this,
				*cluster, temperature);
		std::cout << "Combining Partial Derivative = " << partialDerivatives[index] << std::endl;
	}

	// Load up everything from the dissociating reactants
	size = dissociationConnectivity.size();
	auto reactants = network->getAll();
	for (int i = 0; i < size; i++) {
		// Figure out if this cluster dissociates to another
		if (dissociationConnectivity[i] == 1) {
			cluster = std::dynamic_pointer_cast<PSICluster>(reactants->at(i));
			// Compute the contribution from the cluster
			index = network->getReactantId(*cluster) - 1;
			partialDerivatives[index] += calculateDissociationConstant(*this,
					*cluster, temperature);
		}
	}

	return partialDerivatives;
}

void PSICluster::combineClusters(
		std::shared_ptr<std::vector<std::shared_ptr<Reactant>>>reactants,
		int maxSize, std::string compoundName) {

			std::map<std::string,int> myComposition = getComposition(), secondComposition;
			int numHe, numV, numI, secondNumHe, secondNumV, secondNumI;
			int otherIndex, productIndex;
			std::vector<int> compositionSizes {0,0,0};
			std::shared_ptr<PSICluster> secondCluster, productCluster;
			// Setup the composition variables for this cluster
			numHe = myComposition["He"];
			numV = myComposition["V"];
			numI = myComposition["I"];

			int reactantVecSize = reactants->size();
			for (int i = 0; i < reactantVecSize; i++) {
				// Get the second reactant, its composition and its index
				secondCluster = std::dynamic_pointer_cast <PSICluster> (reactants->at(i));
				secondComposition = secondCluster->getComposition();
				secondNumHe = secondComposition["He"];
				secondNumV = secondComposition["V"];
				secondNumI = secondComposition["I"];
				otherIndex = network->getReactantId(*secondCluster) - 1;
				int productSize = size + secondCluster->getSize();
				// Get and handle product for compounds
				if (compoundName == "HeV" || compoundName == "HeI") {
					// Modify the composition vector
					compositionSizes[0] = numHe + secondNumHe;
					compositionSizes[1] = numV + secondNumV;
					compositionSizes[2] = numI + secondNumI;
					// Get the product
					productCluster = std::dynamic_pointer_cast < PSICluster
					> (network->getCompound(compoundName, compositionSizes));
				} else {
					// Just get the product if it is a single-species
					productCluster = std::dynamic_pointer_cast < PSICluster
					> (network->get(compoundName, productSize));
				}
				// React if the size of the product is valid and it exists in the network
				if (productSize <= maxSize && productCluster) {
					// Setup the connectivity array for the second reactant
					reactionConnectivity[otherIndex] = 1;
					// FIXME! - Debug output
//			std::cout << secondCluster->getSize() << secondCluster->getName() << ": "
//			<< "reactionConnectivity["<< otherIndex << "] = " << reactionConnectivity[otherIndex] << std::endl;
					// Setup the connectivity array for the product
					productIndex = network->getReactantId(*productCluster) - 1;
					reactionConnectivity[productIndex] = 1;
					// FIXME! - Debug output
//			std::cout << productSize << compoundName << ": " << "reactionConnectivity["<< productIndex << "] = "
//			<< reactionConnectivity[productIndex] << std::endl;
					// Push the product onto the list of clusters that combine with this one
					combiningReactants.push_back(secondCluster);
				}
			}

			return;
		}

void PSICluster::replaceInCompound(
		std::shared_ptr<std::vector<std::shared_ptr<Reactant>>>reactants,
		std::string oldComponentName, std::string newComponentName) {

			// Local Declarations
			std::map<std::string, int> myComp = getComposition(),
			secondReactantComp, productReactantComp;
			int myComponentNumber = myComp[newComponentName];
			int productNumHe = 0, productNumV = 0, productNumI = 0;
			int numReactants = reactants->size();
			int secondIndex = 0, productIndex = 0;
			std::shared_ptr<Reactant> secondReactant, productReactant;
			std::vector<int> productCompositionVector(myComp.size());

			// Loop over all of the extra reactants in this reaction and handle the replacement
			for (int i = 0; i < numReactants; i++) {
				// Get the second reactant and its composition
				secondReactant = reactants->at(i);
				secondReactantComp = secondReactant->getComposition();
				// Create the composition vector
				productReactantComp = secondReactantComp;
				// Updated the modified components
				productReactantComp[oldComponentName] =
				secondReactantComp[oldComponentName] - myComponentNumber;
				// Create the composition vector -- FIXME! This should be general!
				productCompositionVector[0] = productReactantComp["He"];
				productCompositionVector[1] = productReactantComp["V"];
				productCompositionVector[2] = productReactantComp["I"];
				// Get the product
				productReactant = network->getCompound(secondReactant->getName(),
						productCompositionVector);
				// If the product exists, mark the proper reaction arrays and it it to the list
				if (productReactant) {
					// Setup the connectivity array for the second reactant
					secondIndex = network->getReactantId(*secondReactant) - 1;
					reactionConnectivity[secondIndex] = 1;
//		std::cout << secondCluster->getSize()
//		<< secondCluster->getName() << ": "
//		<< "reactionConnectivity[" << otherIndex << "] = "
//		<< reactionConnectivity[otherIndex] << std::endl;
					// Setup the connectivity array for the product
					productIndex = network->getReactantId(*productReactant) - 1;
					reactionConnectivity[productIndex] = 1;
					// FIXME! - Debug output
//		std::cout << productSize << compoundName << ": "
//		<< "reactionConnectivity[" << productIndex << "] = "
//		<< reactionConnectivity[productIndex] << std::endl;
					// Push the product onto the list of clusters that combine with this one
					combiningReactants.push_back(secondReactant);
				}
			}

			return;
		}

void PSICluster::fillVWithI(std::string secondClusterName,
		std::shared_ptr<std::vector<std::shared_ptr<Reactant> > > reactants) {
	// Local Declarations
	std::shared_ptr<PSICluster> secondCluster, productCluster;
	std::string firstClusterName = getName(), productClusterName;
	int firstClusterSize = 0, secondClusterSize = 0, productClusterSize = 0;
	int secondIndex = 0, productIndex = 0, reactantVecSize = 0;

	// Get the number of V or I in this cluster (the "first")
	firstClusterSize = (getComposition())[firstClusterName];
	// Look at all of the second clusters, either V or I, and determine
	// if a connection exists.
	reactantVecSize = reactants->size();
	for (int i = 0; i < reactantVecSize; i++) {
		// Get the second cluster its size
		secondCluster = std::dynamic_pointer_cast<PSICluster>(reactants->at(i));
		secondClusterSize =
				(secondCluster->getComposition())[secondClusterName];
		// We have to switch on cluster type to make sure that the annihilation
		// is computed correctly.
		if (firstClusterName == "V") {
			// Compute the product size and set the product name for the case
			// where I is the second cluster
			if (secondClusterSize > firstClusterSize) {
				productClusterSize = secondClusterSize - firstClusterSize;
				productClusterName = "I";
			} else if (secondClusterSize < firstClusterSize) {
				productClusterSize = firstClusterSize - secondClusterSize;
				productClusterName = "V";
			}
		} else if (firstClusterName == "I") {
			// Compute the product size and set the product name for the case
			// where V is the second cluster
			if (firstClusterSize > secondClusterSize) {
				productClusterSize = firstClusterSize - secondClusterSize;
				productClusterName = "I";
			} else if (firstClusterSize < secondClusterSize) {
				productClusterSize = secondClusterSize - firstClusterSize;
				productClusterName = "V";
			}
		}
		// Get the product
		productCluster = std::dynamic_pointer_cast<PSICluster>(
				network->get(productClusterName, productClusterSize));
		// Only deal with this reaction if the product exists. Otherwise the whole reaction is forbidden.
		if (productCluster) {
			// Setup the connectivity array for the second reactant
			secondIndex = network->getReactantId(*secondCluster) - 1;
			reactionConnectivity[secondIndex] = 1;
			// Push the product onto the list of clusters that combine with this one
			combiningReactants.push_back(secondCluster);
			// FIXME! - Debug output
//			std::cout << secondCluster->getSize() << secondCluster->getName() << ": "
//			<< "reactionConnectivity["<< otherIndex << "] = " << reactionConnectivity[otherIndex] << std::endl;
			// Setup the connectivity array for the product
			productIndex = network->getReactantId(*productCluster) - 1;
			reactionConnectivity[productIndex] = 1;
			// FIXME! - Debug output
//			std::cout << productSize << compoundName << ": " << "reactionConnectivity["<< productIndex << "] = "
//			<< reactionConnectivity[productIndex] << std::endl;
		}
	}
}
