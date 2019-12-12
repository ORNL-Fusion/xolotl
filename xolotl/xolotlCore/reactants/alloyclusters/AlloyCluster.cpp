#include <algorithm>
#include "AlloyCluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>
#include <AlloyClusterReactionNetwork.h>

using namespace xolotlCore;

void AlloyCluster::resultFrom(ProductionReaction& reaction,
		IReactant& product) {
	// Cast the reacting clusters
	auto& cluster1 = static_cast<AlloyCluster&>(reaction.first);
	auto& cluster2 = static_cast<AlloyCluster&>(reaction.second);
	auto& prodCluster = static_cast<AlloyCluster&>(product);

	// Compute the overlap
	auto & alloyNetwork = static_cast<AlloyClusterReactionNetwork&>(network);
	int width1 = cluster1.getSectionWidth();
	int size1 = cluster1.getSize();
	int width2 = cluster2.getSectionWidth();
	int size2 = cluster2.getSize();
	int prodWidth = prodCluster.getSectionWidth(), prodSize =
			prodCluster.getSize();
	int lo1 = alloyNetwork.typeSwitch(cluster1.getType())
			* alloyNetwork.typeSwitch(prodCluster.getType())
			* ((int) ((double) size1 - (double) width1 / 2.0) + 1), lo2 =
			alloyNetwork.typeSwitch(cluster2.getType())
					* alloyNetwork.typeSwitch(prodCluster.getType())
					* ((int) ((double) size2 - (double) width2 / 2.0) + 1),
			hi1 = alloyNetwork.typeSwitch(cluster1.getType())
					* alloyNetwork.typeSwitch(prodCluster.getType())
					* ((int) ((double) size1 + (double) width1 / 2.0)), hi2 =
					alloyNetwork.typeSwitch(cluster2.getType())
							* alloyNetwork.typeSwitch(prodCluster.getType())
							* ((int) ((double) size2 + (double) width2 / 2.0));
	int prodLo = ((int) ((double) prodSize - (double) prodWidth / 2.0) + 1),
			prodHi = ((int) ((double) prodSize + (double) prodWidth / 2.0));

	int overlap = std::min(prodHi, hi1 + hi2) - std::max(prodLo, lo1 + lo2) + 1;

	// Skip if the reaction doesn't overlap
	if (overlap < 1)
		return;

	// Create the pair
	ClusterPair pair(reaction, &cluster1, &cluster2);
	// Compute the coefficients
	pair.a0 = overlap;
	if (width1 > 1) {
		pair.a2 = 2.0
				* firstOrderSum(std::max(prodLo - lo2, lo1),
						std::min(prodHi - hi2, hi1), (double) (lo1 + hi1) / 2.0)
				/ (double) (hi1 - lo1);
	}
	if (width2 > 1) {
		pair.a1 = 2.0
				* firstOrderSum(std::max(prodLo - lo1, lo2),
						std::min(prodHi - hi1, hi2), (double) (lo2 + hi2) / 2.0)
				/ (double) ((hi2 - lo2));

	}
	if (width1 > 1 && width2 > 1) {
		// Should never happen for now
		std::cout << "Both reactants are super: " << cluster1.getName() << " + "
				<< cluster2.getName() << " -> " << name << std::endl;
	}

	// Add the pair
	reactingPairs.emplace_back(pair);

	return;
}

void AlloyCluster::resultFrom(ProductionReaction& reaction, int[4], int[4]) {

	// Add a cluster pair for given reaction
	reactingPairs.emplace_back(
			reaction,  // TODO verify this is correct
			&static_cast<AlloyCluster&>(reaction.first),
			&static_cast<AlloyCluster&>(reaction.second));
	// Setup the connectivity array
	setReactionConnectivity(reaction.first.getId());
	setReactionConnectivity(reaction.second.getId());

	return;
}

void AlloyCluster::resultFrom(ProductionReaction& reaction, double *coef) {

	// Add a cluster pair for given reaction
	reactingPairs.emplace_back(
			reaction,  // TODO verify this is correct
			&static_cast<AlloyCluster&>(reaction.first),
			&static_cast<AlloyCluster&>(reaction.second));
	auto& newPair = reactingPairs.back();

	// Setup the connectivity array
	setReactionConnectivity(reaction.first.getId());
	setReactionConnectivity(reaction.second.getId());

	// Update the coefs
	newPair.a0 = coef[0];
	newPair.a1 = coef[1];
	newPair.a2 = coef[2];

	return;
}

void AlloyCluster::participateIn(ProductionReaction& reaction,
		IReactant& product) {
	// Look for the other cluster
	auto& otherCluster = static_cast<AlloyCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);
	auto& prodCluster = static_cast<AlloyCluster&>(product);

	// Compute the overlap
	auto & alloyNetwork = static_cast<AlloyClusterReactionNetwork&>(network);
	int width1 = otherCluster.getSectionWidth();
	int size1 = otherCluster.getSize();
	int width2 = getSectionWidth();
	int size2 = getSize();
	int prodWidth = prodCluster.getSectionWidth(), prodSize =
			prodCluster.getSize();
	int lo1 = alloyNetwork.typeSwitch(otherCluster.getType())
			* alloyNetwork.typeSwitch(prodCluster.getType())
			* ((int) ((double) size1 - (double) width1 / 2.0) + 1), lo2 =
			alloyNetwork.typeSwitch(getType())
					* alloyNetwork.typeSwitch(prodCluster.getType())
					* ((int) ((double) size2 - (double) width2 / 2.0) + 1),
			hi1 = alloyNetwork.typeSwitch(otherCluster.getType())
					* alloyNetwork.typeSwitch(prodCluster.getType())
					* ((int) ((double) size1 + (double) width1 / 2.0)), hi2 =
					alloyNetwork.typeSwitch(getType())
							* alloyNetwork.typeSwitch(prodCluster.getType())
							* ((int) ((double) size2 + (double) width2 / 2.0));
	int prodLo = ((int) ((double) prodSize - (double) prodWidth / 2.0) + 1),
			prodHi = ((int) ((double) prodSize + (double) prodWidth / 2.0));

	int overlap = std::min(prodHi, hi1 + hi2) - std::max(prodLo, lo1 + lo2) + 1;

	// Skip if the reaction doesn't overlap
	if (overlap < 1)
		return;

	// Create the pair
	CombiningCluster pair(reaction, &otherCluster);

	// Compute the coefficients
	pair.a0 = overlap;

	if (width1 > 1) {
		pair.a1 = 2.0
				* firstOrderSum(std::max(prodLo - lo2, lo1),
						std::min(prodHi - hi2, hi1), (double) (lo1 + hi1) / 2.0)
				/ (double) (hi1 - lo1);
	}

	// Add the combining cluster to list of clusters that combine with us
	combiningReactants.emplace_back(pair);

	return;
}

void AlloyCluster::participateIn(ProductionReaction& reaction, int[4]) {
	// Look for the other cluster
	auto& otherCluster = static_cast<AlloyCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);


	// Create the pair
	CombiningCluster pair(reaction, &otherCluster);
	// Compute the coefficients
	// This is called by the clusters combining together
	pair.a0 = 1.0;
	// Add the combining cluster to list of clusters that combine with us
	combiningReactants.emplace_back(pair);

	// Setup the connectivity array
	setReactionConnectivity(id);
	setReactionConnectivity(otherCluster.getId());

	return;
}

void AlloyCluster::participateIn(ProductionReaction& reaction, double *coef) {
	// Look for the other cluster
	auto& otherCluster = static_cast<AlloyCluster&>(
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

	// Update the coefs
	newComb.a0 = coef[0];
	newComb.a1 = coef[1];

	return;
}

void AlloyCluster::participateIn(DissociationReaction& reaction,
		IReactant& disso) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<AlloyCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);
	auto& dissoCluster = static_cast<AlloyCluster&>(disso);

	// Compute the overlap
	auto & alloyNetwork = static_cast<AlloyClusterReactionNetwork&>(network);
	int width1 = emittedCluster.getSectionWidth();
	int size1 = emittedCluster.getSize();
	int width2 = getSectionWidth();
	int size2 = getSize();
	int dissoWidth = dissoCluster.getSectionWidth(), dissoSize =
			dissoCluster.getSize();
	int lo1 = alloyNetwork.typeSwitch(emittedCluster.getType())
			* alloyNetwork.typeSwitch(dissoCluster.getType())
			* ((int) ((double) size1 - (double) width1 / 2.0) + 1), lo2 =
			alloyNetwork.typeSwitch(getType())
					* alloyNetwork.typeSwitch(dissoCluster.getType())
					* ((int) ((double) size2 - (double) width2 / 2.0) + 1),
			hi1 = alloyNetwork.typeSwitch(emittedCluster.getType())
					* alloyNetwork.typeSwitch(dissoCluster.getType())
					* ((int) ((double) size1 + (double) width1 / 2.0)), hi2 =
					alloyNetwork.typeSwitch(getType())
							* alloyNetwork.typeSwitch(dissoCluster.getType())
							* ((int) ((double) size2 + (double) width2 / 2.0));
	int dissoLo = ((int) ((double) dissoSize - (double) dissoWidth / 2.0) + 1),
			dissoHi = ((int) ((double) dissoSize + (double) dissoWidth / 2.0));

	int overlap = std::min(dissoHi, hi1 + hi2) - std::max(dissoLo, lo1 + lo2)
			+ 1;

	// Skip if the reaction doesn't overlap
	if (overlap < 1)
		return;

	// Create the pair
	ClusterPair pair(reaction, &dissoCluster, &emittedCluster);
	// Compute the coefficients
	pair.a0 = overlap;
	if (dissoWidth > 1) {
		pair.a1 = 2.0
				* firstOrderSum(std::max(dissoLo, lo1 + lo2),
						std::min(dissoHi, hi1 + hi2),
						(double) (dissoLo + dissoHi) / 2.0)
				/ (double) (dissoHi - dissoLo);
	}

	// Add the pair to the vector
	dissociatingPairs.emplace_back(pair);

	return;
}

void AlloyCluster::participateIn(DissociationReaction& reaction, int[4],
		int[4]) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<AlloyCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Add a pair where it is important that the
	// dissociating cluster is the first one
	dissociatingPairs.emplace_back(
			reaction,  // TODO is this correct?
			&static_cast<AlloyCluster&>(reaction.dissociating),
			&static_cast<AlloyCluster&>(emittedCluster));

	// Setup the connectivity array
	setDissociationConnectivity(reaction.dissociating.getId());

	return;
}

void AlloyCluster::participateIn(DissociationReaction& reaction, double *coef) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<AlloyCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Add a pair where it is important that the
	// dissociating cluster is the first one
	dissociatingPairs.emplace_back(
			reaction,  // TODO is this correct?
			&static_cast<AlloyCluster&>(reaction.dissociating),
			&static_cast<AlloyCluster&>(emittedCluster));
	auto& newPair = dissociatingPairs.back();

	// Setup the connectivity array
	setDissociationConnectivity(reaction.dissociating.getId());

	// Set the coefs
	newPair.a0 = coef[0];
	newPair.a1 = coef[1];

	return;
}

void AlloyCluster::emitFrom(DissociationReaction& reaction, IReactant& disso) {
	// Cast the reacting clusters
	auto& cluster1 = static_cast<AlloyCluster&>(reaction.first);
	auto& cluster2 = static_cast<AlloyCluster&>(reaction.second);
	auto& dissoCluster = static_cast<AlloyCluster&>(disso);

	// Compute the overlap
	auto & alloyNetwork = static_cast<AlloyClusterReactionNetwork&>(network);
	int width1 = cluster1.getSectionWidth();
	int size1 = cluster1.getSize();
	int width2 = cluster2.getSectionWidth();
	int size2 = cluster2.getSize();
	int dissoWidth = dissoCluster.getSectionWidth(), dissoSize =
			dissoCluster.getSize();
	int lo1 = alloyNetwork.typeSwitch(cluster1.getType())
			* alloyNetwork.typeSwitch(dissoCluster.getType())
			* ((int) ((double) size1 - (double) width1 / 2.0) + 1), lo2 =
			alloyNetwork.typeSwitch(cluster2.getType())
					* alloyNetwork.typeSwitch(dissoCluster.getType())
					* ((int) ((double) size2 - (double) width2 / 2.0) + 1),
			hi1 = alloyNetwork.typeSwitch(cluster1.getType())
					* alloyNetwork.typeSwitch(dissoCluster.getType())
					* ((int) ((double) size1 + (double) width1 / 2.0)), hi2 =
					alloyNetwork.typeSwitch(cluster2.getType())
							* alloyNetwork.typeSwitch(dissoCluster.getType())
							* ((int) ((double) size2 + (double) width2 / 2.0));
	int dissoLo = ((int) ((double) dissoSize - (double) dissoWidth / 2.0) + 1),
			dissoHi = ((int) ((double) dissoSize + (double) dissoWidth / 2.0));

	int overlap = std::min(dissoHi, hi1 + hi2) - std::max(dissoLo, lo1 + lo2)
			+ 1;

	// Skip if the reaction doesn't overlap
	if (overlap < 1)
		return;

	// Create the pair
	ClusterPair pair(reaction, &cluster1, &cluster2);

	// Compute the coefficients
	pair.a0 = overlap;

	// Add the pair
	emissionPairs.emplace_back(pair);

	return;
}

void AlloyCluster::emitFrom(DissociationReaction& reaction, int[4]) {

	// Add the pair of emitted clusters.
	emissionPairs.emplace_back(
			reaction, // TODO is this correct?
			&static_cast<AlloyCluster&>(reaction.first),
			&static_cast<AlloyCluster&>(reaction.second));

	// Setup the connectivity array to itself
	setReactionConnectivity(id);

	return;
}

void AlloyCluster::emitFrom(DissociationReaction& reaction, double *coef) {

	// Add the pair of emitted clusters.
	emissionPairs.emplace_back(
			reaction, // TODO is this correct?
			&static_cast<AlloyCluster&>(reaction.first),
			&static_cast<AlloyCluster&>(reaction.second));
	auto& newPair = emissionPairs.back();

	// Setup the connectivity array to itself
	setReactionConnectivity(id);

	// Set the coefs
	newPair.a0 = coef[0];

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

std::vector<int> AlloyCluster::getReactionConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(reactionConnectivitySet, network.getDOF());
}

std::vector<int> AlloyCluster::getDissociationConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(dissociationConnectivitySet,
			network.getDOF());
}

void AlloyCluster::resetConnectivities() {
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
				AlloyCluster const& combCluster = *cc.combining;
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

double AlloyCluster::getMoment() const {
	return 0.0;
}

double AlloyCluster::getTotalFlux(int i) {
	// Get the fluxes
	double prodFlux = getProductionFlux(i);
	double dissFlux = getDissociationFlux(i);
	double combFlux = getCombinationFlux(i);
	double emissFlux = getEmissionFlux(i);

	return prodFlux - combFlux + dissFlux - emissFlux;
}

double AlloyCluster::getDissociationFlux(int xi) const {

	// Sum dissociation flux over all pairs that dissociate to form this one.
	double flux =
			std::accumulate(dissociatingPairs.begin(), dissociatingPairs.end(),
					0.0, [&xi](double running, const ClusterPair& currPair) {
						// Get the dissociating clusters
					AlloyCluster* dissociatingCluster = currPair.first;
					double l0A = dissociatingCluster->getConcentration();
					double l1A = dissociatingCluster->getMoment();
					// Update the flux
					return running + (currPair.reaction.kConstant[xi] * (currPair.a0 * l0A + currPair.a1 * l1A));
				});

	// Return the flux
	return flux;
}

double AlloyCluster::getEmissionFlux(int xi) const {

	// Sum reaction rate constants over all emission pair reactions.
	double flux = std::accumulate(emissionPairs.begin(), emissionPairs.end(),
			0.0, [&xi](double running, const ClusterPair& currPair) {
				return running + currPair.reaction.kConstant[xi] * currPair.a0;
			});

	return flux * concentration;
}

double AlloyCluster::getProductionFlux(int xi) const {
	// Local declarations
	double flux = 0.0;

	// Sum over all the reacting pairs
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[&flux,&xi](ClusterPair const& currPair) {
				// Get the two reacting clusters
				AlloyCluster* firstReactant = currPair.first;
				AlloyCluster* secondReactant = currPair.second;
				// We know the first one is always the single one
				double l0A = firstReactant->getConcentration();
				double l1A = firstReactant->getMoment();
				double l0B = secondReactant->getConcentration();
				double l1B = secondReactant->getMoment();
				// Update the flux
				// The double moment term is not possible because 2 super can't react together
				flux += currPair.reaction.kConstant[xi] *
						(currPair.a0 * l0A * l0B + currPair.a1 * l0B * l1A + currPair.a2 * l0A * l1B);
			});

	// Return the production flux
	return flux;
}

double AlloyCluster::getCombinationFlux(int xi) const {

	double flux = std::accumulate(combiningReactants.begin(),
			combiningReactants.end(), 0.0,
			[xi, this](double running, const CombiningCluster& currPair) {
				// Get the cluster that combines with this one
				AlloyCluster const& combiningCluster = *currPair.combining;
				double l0A = combiningCluster.getConcentration();
				double l1A = combiningCluster.getMoment();
				// Update the flux
				return running +
				(currPair.reaction.kConstant[xi] *
						(currPair.a0 * l0A + currPair.a1 * l1A));
			});

	return flux * concentration;
}

std::vector<double> AlloyCluster::getPartialDerivatives(int i) const {
	// Local Declarations
	std::vector<double> partials(network.getDOF(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials, i);
	getCombinationPartialDerivatives(partials, i);
	getDissociationPartialDerivatives(partials, i);
	getEmissionPartialDerivatives(partials, i);

	return partials;
}

void AlloyCluster::getPartialDerivatives(std::vector<double> & partials,
		int i) const {
	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials, i);
	getCombinationPartialDerivatives(partials, i);
	getDissociationPartialDerivatives(partials, i);
	getEmissionPartialDerivatives(partials, i);

	return;
}

void AlloyCluster::getProductionPartialDerivatives(
		std::vector<double> & partials, int xi) const {

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[&partials,&xi](ClusterPair const& currPair) {
				// Get the two reacting clusters
				AlloyCluster* firstReactant = currPair.first;
				AlloyCluster* secondReactant = currPair.second;
				double l0A = firstReactant->getConcentration();
				double l1A = firstReactant->getMoment();
				double l0B = secondReactant->getConcentration();
				double l1B = secondReactant->getMoment();

				// Compute the contribution from the first part of the reacting pair
				auto value = currPair.reaction.kConstant[xi];
				auto index = firstReactant->id - 1;
				partials[index] += value * (currPair.a0 * l0B + currPair.a2 * l1B);
				index = firstReactant->momId[0] - 1;
				partials[index] += value * currPair.a1 * l0B;
				// Compute the contribution from the second part of the reacting pair
				index = secondReactant->id - 1;
				partials[index] += value * (currPair.a0 * l0A + currPair.a1 * l1A);
				index = secondReactant->momId[0] - 1;
				partials[index] += value * currPair.a2 * l0A;
			});

	return;
}

void AlloyCluster::getCombinationPartialDerivatives(
		std::vector<double> & partials, int xi) const {

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[this,&partials,&xi](const CombiningCluster& cc) {
				AlloyCluster const& cluster = *cc.combining;
				Reaction const& currReaction = cc.reaction;
				double l0A = cluster.getConcentration();
				double l1A = cluster.getMoment();

				// Remember that the flux due to combinations is OUTGOING (-=)!
				// Compute the contribution from this cluster
				partials[id - 1] -= currReaction.kConstant[xi] * (cc.a0 * l0A + cc.a1 * l1A);
				// Compute the contribution from the combining cluster
				double value = currReaction.kConstant[xi] * concentration;

				partials[cluster.id - 1] -= value * cc.a0;
				partials[cluster.momId[0] - 1] -= value * cc.a1;
			});

	return;
}

void AlloyCluster::getDissociationPartialDerivatives(
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
				AlloyCluster* cluster = currPair.first;
				Reaction const& currReaction = currPair.reaction;
				partials[cluster->id - 1] += currReaction.kConstant[xi] * currPair.a0;
				partials[cluster->momId[0] - 1] += currReaction.kConstant[xi] *
				currPair.a1;
			});

	return;
}

void AlloyCluster::getEmissionPartialDerivatives(std::vector<double> & partials,
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
				return running + currPair.reaction.kConstant[xi] * currPair.a0;
			});

	// Recall emission flux is OUTGOING
	partials[id - 1] -= emissionFlux;

	return;
}

double AlloyCluster::getLeftSideRate(int i) const {

	// Sum reaction rate contributions over all combining clusters.
	double combiningRateTotal = std::accumulate(combiningReactants.begin(),
			combiningReactants.end(), 0.0,
			[&i](double running, const CombiningCluster& currPair) {
				AlloyCluster const& cluster = *currPair.combining;
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

std::vector<std::vector<double> > AlloyCluster::getProdVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the reacting pairs
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[&toReturn](const ClusterPair& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.second->getId() - 1);
				tempVec.push_back(currPair.a0);
				tempVec.push_back(currPair.a1);
				tempVec.push_back(currPair.a2);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > AlloyCluster::getCombVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the combining reactants
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[&toReturn](const CombiningCluster& cc) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(cc.combining->getId() - 1);
				tempVec.push_back(cc.a0);
				tempVec.push_back(cc.a1);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > AlloyCluster::getDissoVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the dissociating pairs
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[&toReturn](const ClusterPair& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.second->getId() - 1);
				tempVec.push_back(currPair.a0);
				tempVec.push_back(currPair.a1);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > AlloyCluster::getEmitVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the emitting pairs
	std::for_each(emissionPairs.begin(), emissionPairs.end(),
			[&toReturn](const ClusterPair& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.second->getId() - 1);
				tempVec.push_back(currPair.a0);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<int> AlloyCluster::getConnectivity() const {
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
