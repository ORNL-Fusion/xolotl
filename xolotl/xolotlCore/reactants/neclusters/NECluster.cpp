#include <algorithm>
#include <cassert>
#include "NECluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

void NECluster::resultFrom(ProductionReaction& reaction, int[4], int[4]) {

	// Add a cluster pair for given reaction 
	reactingPairs.emplace_back(
			reaction,  // TODO verify this is correct
			&static_cast<NECluster&>(reaction.first),
			&static_cast<NECluster&>(reaction.second));
	// Setup the connectivity array
	setReactionConnectivity(reaction.first.getId());
	setReactionConnectivity(reaction.second.getId());

	return;
}

void NECluster::resultFrom(ProductionReaction& reaction, double *coef) {

	// Add a cluster pair for given reaction
	reactingPairs.emplace_back(
			reaction,  // TODO verify this is correct
			&static_cast<NECluster&>(reaction.first),
			&static_cast<NECluster&>(reaction.second));
	auto& newPair = reactingPairs.back();

	// Setup the connectivity array
	setReactionConnectivity(reaction.first.getId());
	setReactionConnectivity(reaction.second.getId());

	// Update the distances
	newPair.firstDistance = coef[0];
	newPair.secondDistance = coef[1];

	return;
}

void NECluster::participateIn(ProductionReaction& reaction, int[4]) {
	// Look for the other cluster
	auto& otherCluster = static_cast<NECluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);
//	// Build a production reaction for it
//	std::unique_ptr<ProductionReaction> newReaction(
//			new ProductionReaction(otherCluster, *this));
//	// Add it to the network
//	auto& prref = network.add(std::move(newReaction));

	// Add the combining cluster to list of clusters that combine with us
	combiningReactants.emplace_back(reaction, &otherCluster);

	// Setup the connectivity array
	setReactionConnectivity(id);
	setReactionConnectivity(otherCluster.getId());

	return;
}

void NECluster::participateIn(ProductionReaction& reaction, double *coef) {
	// Look for the other cluster
	auto& otherCluster = static_cast<NECluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);
//	// Build a production reaction for it
//	std::unique_ptr<ProductionReaction> newReaction(
//			new ProductionReaction(otherCluster, *this));
//	// Add it to the network
//	auto& prref = network.add(std::move(newReaction));

	// Add the combining cluster to list of clusters that combine with us
	combiningReactants.emplace_back(reaction, &otherCluster);
	auto& newComb = combiningReactants.back();

	// Setup the connectivity array
	setReactionConnectivity(id);
	setReactionConnectivity(otherCluster.getId());

	// Update the distances
	newComb.distance = coef[0];

	return;
}

void NECluster::participateIn(DissociationReaction& reaction, int[4], int[4]) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<NECluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Add a pair where it is important that the
	// dissociating cluster is the first one
	dissociatingPairs.emplace_back(
			reaction,  // TODO is this correct?
			&static_cast<NECluster&>(reaction.dissociating),
			&static_cast<NECluster&>(emittedCluster));

	// Setup the connectivity array
	setDissociationConnectivity(reaction.dissociating.getId());

	return;
}

void NECluster::participateIn(DissociationReaction& reaction, double *coef) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<NECluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Add a pair where it is important that the
	// dissociating cluster is the first one
	dissociatingPairs.emplace_back(
			reaction,  // TODO is this correct?
			&static_cast<NECluster&>(reaction.dissociating),
			&static_cast<NECluster&>(emittedCluster));
	auto& newPair = dissociatingPairs.back();

	// Setup the connectivity array
	setDissociationConnectivity(reaction.dissociating.getId());

	// Set the distance
	newPair.firstDistance = coef[0];

	return;
}

void NECluster::emitFrom(DissociationReaction& reaction, int[4]) {

	// Add the pair of emitted clusters.
	emissionPairs.emplace_back(
			reaction, // TODO is this correct?
			&static_cast<NECluster&>(reaction.first),
			&static_cast<NECluster&>(reaction.second));

	// Setup the connectivity array to itself
	setReactionConnectivity(id);

	return;
}

void NECluster::emitFrom(DissociationReaction& reaction, double *coef) {

	// Add the pair of emitted clusters.
	emissionPairs.emplace_back(
			reaction, // TODO is this correct?
			&static_cast<NECluster&>(reaction.first),
			&static_cast<NECluster&>(reaction.second));

	// Setup the connectivity array to itself
	setReactionConnectivity(id);

	// Nothing more to do

	return;
}

void NECluster::optimizeReactions() {
	// Loop on the pairs to add reactions to the network
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[this](ClusterPair& currPair) {
				// Create the corresponding production reaction
				std::unique_ptr<ProductionReaction> newReaction(new ProductionReaction(*currPair.first, *currPair.second));
				// Add it to the network
				auto& prref = network.add(std::move(newReaction));
				// Link it to the pair
				currPair.reaction = prref;
			});

	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[this](CombiningCluster& cc) {
				// Create the corresponding production reaction
				std::unique_ptr<ProductionReaction> newReaction(new ProductionReaction(*cc.combining, *this));
				// Add it to the network
				auto& prref = network.add(std::move(newReaction));
				// Link it to the pair
				cc.reaction = prref;
			});

	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[this](ClusterPair& currPair) {
				// Create the corresponding dissociation reaction
				std::unique_ptr<DissociationReaction> newReaction(new DissociationReaction(*currPair.first, *currPair.second, *this));
				// Add it to the network
				auto& drref = network.add(std::move(newReaction));
				// Create the corresponding reverse reaction
				std::unique_ptr<ProductionReaction> newReverseReaction(new ProductionReaction(*currPair.second, *this));
				// Add it to the network
				auto& prref = network.add(std::move(newReverseReaction));
				// Link it
				drref.reverseReaction = &prref;
				// Link it to the pair
				currPair.reaction = drref;
			});

	std::for_each(emissionPairs.begin(), emissionPairs.end(),
			[this](ClusterPair& currPair) {
				// Create the corresponding dissociation reaction
				std::unique_ptr<DissociationReaction> newReaction(new DissociationReaction(*this, *currPair.first, *currPair.second));
				// Add it to the network
				auto& drref = network.add(std::move(newReaction));
				// Create the corresponding reverse reaction
				std::unique_ptr<ProductionReaction> newReverseReaction(new ProductionReaction(*currPair.first, *currPair.second));
				// Add it to the network
				auto& prref = network.add(std::move(newReverseReaction));
				// Link it
				drref.reverseReaction = &prref;
				// Link it to the pair
				currPair.reaction = drref;
			});

	return;
}

static std::vector<int> getFullConnectivityVector(std::set<int> connectivitySet,
		int size) {
	// Create a vector of zeroes with size equal to the network size
	std::vector<int> connectivity(size);

	// Set the value of the connectivity array to one for each element that is
	// in the set.
	for (auto it = connectivitySet.begin(); it != connectivitySet.end(); ++it) {
		connectivity[(*it) - 1] = 1;
	}

	return connectivity;
}

std::vector<int> NECluster::getReactionConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(reactionConnectivitySet, network.getDOF());
}

std::vector<int> NECluster::getDissociationConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(dissociationConnectivitySet,
			network.getDOF());
}

void NECluster::resetConnectivities() {
	// Shrink the arrays to save some space. (About 10% or so.)
	reactingPairs.shrink_to_fit();
	combiningReactants.shrink_to_fit();
	dissociatingPairs.shrink_to_fit();
	emissionPairs.shrink_to_fit();

	// Clear both sets
	reactionConnectivitySet.clear();
	dissociationConnectivitySet.clear();

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);
	setDissociationConnectivity(id);
	setReactionConnectivity(momId[0]);
	setDissociationConnectivity(momId[0]);

	// Apply to each reacting pair.
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[this](const ClusterPair& currPair) {
				// The cluster is connecting to both clusters in the pair
				setReactionConnectivity(currPair.first->id);
				setReactionConnectivity(currPair.first->momId[0]);
				setReactionConnectivity(currPair.second->id);
				setReactionConnectivity(currPair.second->momId[0]);
			});

	// Apply to each combining cluster.
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[this](const CombiningCluster& cc) {
				// The cluster is connecting to the combining cluster
				NECluster const& combCluster = *cc.combining;
				setReactionConnectivity(combCluster.id);
				setReactionConnectivity(combCluster.momId[0]);
			});

	// Apply to each effective dissociating pair
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[this](const ClusterPair& currPair) {
				// The cluster is connecting to the dissociating cluster which
				// is the first one by definition
				setDissociationConnectivity(currPair.first->id);
				setDissociationConnectivity(currPair.first->momId[0]);
			});

	// Don't apply to the emission pairs because
	// this cluster is not connected to them

	return;
}

void NECluster::updateFromNetwork() {

	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	return;
}

double NECluster::getMoment() const {
	return 0.0;
}

double NECluster::getTotalFlux(int i) {
	// Get the fluxes
	double prodFlux = getProductionFlux(i);
	double dissFlux = getDissociationFlux(i);
	double combFlux = getCombinationFlux(i);
	double emissFlux = getEmissionFlux(i);

	return prodFlux - combFlux + dissFlux - emissFlux;
}

double NECluster::getDissociationFlux(int xi) {

	// Sum dissociation flux over all pairs that dissociate to form this one.
	double flux =
			std::accumulate(dissociatingPairs.begin(), dissociatingPairs.end(),
					0.0, [&xi](double running, const ClusterPair& currPair) {
						// Get the dissociating cluster
					auto& dissociatingCluster = currPair.first;
					// Calculate the Dissociation flux
					Reaction const& currReaction = currPair.reaction;
					return running + (currReaction.kConstant[xi] *
							dissociatingCluster->getConcentration(currPair.firstDistance));
				});

	// Return the flux
	return flux;
}

double NECluster::getEmissionFlux(int xi) {

	// Sum reaction rate constants over all emission pair reactions.
	double flux = std::accumulate(emissionPairs.begin(), emissionPairs.end(),
			0.0, [&xi](double running, const ClusterPair& currPair) {
				Reaction const& currReaction = currPair.reaction;
				return running + currReaction.kConstant[xi];
			});

	return flux * concentration;
}

double NECluster::getProductionFlux(int xi) {
	// Local declarations
	double flux = 0.0;

	// Sum over all the reacting pairs
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[&flux,&xi](ClusterPair const& currPair) {
				// Get the two reacting clusters
				NECluster* firstReactant = currPair.first;
				NECluster* secondReactant = currPair.second;
				// Update the flux
				Reaction const& currReaction = currPair.reaction;
				flux += currReaction.kConstant[xi]
				* firstReactant->getConcentration(
						currPair.firstDistance)
				* secondReactant->getConcentration(
						currPair.secondDistance);
			});

	// Return the production flux
	return flux;
}

double NECluster::getCombinationFlux(int xi) {

	double flux = std::accumulate(combiningReactants.begin(),
			combiningReactants.end(), 0.0,
			[xi](double running, const CombiningCluster& currPair) {
				// Get the cluster that combines with this one
				NECluster const& combiningCluster = *currPair.combining;
				Reaction const& currReaction = currPair.reaction;

				// Calculate Second term of production flux
				return running +
				(currReaction.kConstant[xi] *
						combiningCluster.getConcentration(currPair.distance));

			});

	return flux * concentration;
}

std::vector<double> NECluster::getPartialDerivatives(int i) const {
	// Local Declarations
	std::vector<double> partials(network.getDOF(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials, i);
	getCombinationPartialDerivatives(partials, i);
	getDissociationPartialDerivatives(partials, i);
	getEmissionPartialDerivatives(partials, i);

	return partials;
}

void NECluster::getPartialDerivatives(std::vector<double> & partials,
		int i) const {
	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials, i);
	getCombinationPartialDerivatives(partials, i);
	getDissociationPartialDerivatives(partials, i);
	getEmissionPartialDerivatives(partials, i);

	return;
}

void NECluster::getProductionPartialDerivatives(std::vector<double> & partials,
		int xi) const {

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[&partials,&xi](ClusterPair const& currPair) {

				Reaction const& currReaction = currPair.reaction;

				// Compute the contribution from the first part of the reacting pair
				auto value = currReaction.kConstant[xi]
				* currPair.second->getConcentration(
						currPair.secondDistance);
				auto index = currPair.first->id - 1;
				partials[index] += value;
				index = currPair.first->momId[0] - 1;
				partials[index] += value * currPair.firstDistance;
				// Compute the contribution from the second part of the reacting pair
				value = currReaction.kConstant[xi]
				* currPair.first->getConcentration(
						currPair.firstDistance);
				index = currPair.second->id - 1;
				partials[index] += value;
				index = currPair.second->momId[0] - 1;
				partials[index] += value * currPair.secondDistance;
			});

	return;
}

void NECluster::getCombinationPartialDerivatives(std::vector<double> & partials,
		int xi) const {

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[this,&partials,&xi](const CombiningCluster& cc) {

				NECluster const& cluster = *cc.combining;
				Reaction const& currReaction = cc.reaction;

				// Remember that the flux due to combinations is OUTGOING (-=)!
				// Compute the contribution from this cluster
				partials[id - 1] -= currReaction.kConstant[xi] *
				cluster.getConcentration(cc.distance);
				// Compute the contribution from the combining cluster
				double value = currReaction.kConstant[xi] * concentration;

				partials[cluster.id - 1] -= value;
				partials[cluster.momId[0] - 1] -= value * cc.distance;
			});

	return;
}

void NECluster::getDissociationPartialDerivatives(
		std::vector<double> & partials, int xi) const {

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[&partials,&xi](const ClusterPair& currPair) {
				// Get the dissociating cluster
				NECluster* cluster = currPair.first;
				Reaction const& currReaction = currPair.reaction;
				partials[cluster->id - 1] += currReaction.kConstant[xi];
				partials[cluster->momId[0] - 1] += currReaction.kConstant[xi] *
				currPair.firstDistance;
			});

	return;
}

void NECluster::getEmissionPartialDerivatives(std::vector<double> & partials,
		int xi) const {

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)
	double emissionFlux = std::accumulate(emissionPairs.begin(),
			emissionPairs.end(), 0.0,
			[&xi](double running, const ClusterPair& currPair) {
				Reaction const& currReaction = currPair.reaction;
				return running + currReaction.kConstant[xi];
			});

	// Recall emission flux is OUTGOING
	partials[id - 1] -= emissionFlux;

	return;
}

double NECluster::getLeftSideRate(int i) const {

	// Sum reaction rate contributions over all combining clusters.
	double combiningRateTotal = std::accumulate(combiningReactants.begin(),
			combiningReactants.end(), 0.0,
			[&i](double running, const CombiningCluster& currPair) {
				NECluster const& cluster = *currPair.combining;
				Reaction const& currReaction = currPair.reaction;

				return running + (currReaction.kConstant[i] *
						cluster.concentration);
			});

	// Sum reaction rate constants over all emission pairs.
	double emissionRateTotal = std::accumulate(emissionPairs.begin(),
			emissionPairs.end(), 0.0,
			[&i](double running, const ClusterPair& currPair) {
				Reaction const& currReaction = currPair.reaction;
				return running + currReaction.kConstant[i];
			});

	return combiningRateTotal + emissionRateTotal;
}

std::vector<std::vector<double> > NECluster::getProdVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the reacting pairs
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[&toReturn](const ClusterPair& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.second->getId() - 1);
				tempVec.push_back(currPair.firstDistance);
				tempVec.push_back(currPair.secondDistance);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > NECluster::getCombVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the combining reactants
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[&toReturn](const CombiningCluster& cc) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(cc.combining->getId() - 1);
				tempVec.push_back(cc.distance);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > NECluster::getDissoVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the dissociating pairs
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[&toReturn](const ClusterPair& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.second->getId() - 1);
				tempVec.push_back(currPair.firstDistance);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > NECluster::getEmitVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the emitting pairs
	std::for_each(emissionPairs.begin(), emissionPairs.end(),
			[&toReturn](const ClusterPair& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.second->getId() - 1);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<int> NECluster::getConnectivity() const {
	int connectivityLength = network.getDOF();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);
	auto reactionConnVector = getReactionConnectivity();
	auto dissociationConnVector = getDissociationConnectivity();

	// The reaction and dissociation vectors must have a length equal to the
	// number of clusters
	if (reactionConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The reaction vector has an incorrect length");
	}
	if (dissociationConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The dissociation vector has an incorrect length");
	}

	// Merge the two vectors such that the final vector contains
	// a 1 at a position if either of the connectivity arrays
	// have a 1
	for (int i = 0; i < connectivityLength; i++) {
		// Consider each connectivity array only if its type is enabled
		connectivity[i] = reactionConnVector[i] || dissociationConnVector[i];
	}

	return connectivity;
}
