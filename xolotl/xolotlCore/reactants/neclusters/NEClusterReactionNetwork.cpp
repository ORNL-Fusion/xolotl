#include "NEClusterReactionNetwork.h"
#include "NECluster.h"
#include "NESuperCluster.h"
#include <xolotlPerf.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <Constants.h>

using namespace xolotlCore;

NEClusterReactionNetwork::NEClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork( { ReactantType::V, ReactantType::I, ReactantType::Xe,
				ReactantType::XeV, ReactantType::XeI, ReactantType::NESuper },
				ReactantType::NESuper, registry) {

	// Initialize default properties
	dissociationsEnabled = true;

	return;
}

double NEClusterReactionNetwork::calculateDissociationConstant(
		const DissociationReaction& reaction) const {

	// If the dissociations are not allowed
	if (!dissociationsEnabled)
		return 0.0;

	// Compute the atomic volume
	double atomicVolume = 0.5 * xolotlCore::uraniumDioxydeLatticeConstant
			* xolotlCore::uraniumDioxydeLatticeConstant
			* xolotlCore::uraniumDioxydeLatticeConstant;

	// Get the rate constant from the reverse reaction
	double kPlus = reaction.reverseReaction->kConstant;

	// Calculate and return
	double bindingEnergy = computeBindingEnergy(reaction);
	double k_minus_exp = exp(
			-1.0 * bindingEnergy / (xolotlCore::kBoltzmann * temperature));
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;

	return k_minus;
}

void NEClusterReactionNetwork::createReactionConnectivity() {
	// Initial declarations
	int firstSize = 0, secondSize = 0, productSize = 0;

	// Single species clustering (Xe)
	// We know here that only Xe_1 can cluster so we simplify the search
	// Xe_(a-i) + Xe_i --> Xe_a
	firstSize = 1;
	auto singleXeCluster = get(Species::Xe, firstSize);
	// Get all the Xe clusters
	auto const& xeClusters = getAll(ReactantType::Xe);
	// Consider each Xe cluster.
	for (auto const& currMapItem : getAll(ReactantType::Xe)) {

		auto& xeReactant = *(currMapItem.second);

		// Get the size of the second reactant and product
		secondSize = xeReactant.getSize();
		productSize = firstSize + secondSize;
		// Get the product cluster for the reaction
		auto product = get(Species::Xe, productSize);
		// Check that the reaction can occur
		if (product
				&& (singleXeCluster->getDiffusionFactor() > 0.0
						|| xeReactant.getDiffusionFactor() > 0.0)) {
			// Create a production reaction
			auto reaction = std::make_shared<ProductionReaction>(
					*singleXeCluster, xeReactant);
			// Tell the reactants that they are in this reaction
			singleXeCluster->participateIn(*reaction);
			xeReactant.participateIn(*reaction);
			product->resultFrom(*reaction);

			// Check if the reverse reaction is allowed
			checkForDissociation(product, reaction);
		}
	}

	return;
}

void NEClusterReactionNetwork::checkForDissociation(
		IReactant * emittingReactant,
		std::shared_ptr<ProductionReaction> reaction) {
	// Check if at least one of the potentially emitted cluster is size one
	if (reaction->first.getSize() != 1 && reaction->second.getSize() != 1) {
		// Don't add the reverse reaction
		return;
	}

	// The reaction can occur, create the dissociation
	// Create a dissociation reaction
	// TODO can this be on the stack?
	std::unique_ptr<DissociationReaction> dissociationReaction(
			new DissociationReaction(*emittingReactant, reaction->first,
					reaction->second));
	// Set the reverse reaction
	dissociationReaction->reverseReaction = reaction.get();
	// Tell the reactants that their are in this reaction
	reaction->first.participateIn(*dissociationReaction);
	reaction->second.participateIn(*dissociationReaction);
	emittingReactant->emitFrom(*dissociationReaction);

	return;
}

void NEClusterReactionNetwork::setTemperature(double temp) {
	ReactionNetwork::setTemperature(temp);

	computeRateConstants();

	return;
}

void NEClusterReactionNetwork::reinitializeNetwork() {
	// Reset the Ids
	int id = 0;
	std::for_each(allReactants.begin(), allReactants.end(),
			[&id](IReactant& currReactant) {
				id++;
				currReactant.setId(id);
				currReactant.setMomentId(id);

				currReactant.optimizeReactions();
			});

	// Get all the super clusters and loop on them
	for (auto const& currMapItem : clusterTypeMap.at(ReactantType::NESuper)) {

		auto& currCluster = static_cast<NESuperCluster&>(*(currMapItem.second));
		id++;
		currCluster.setMomentId(id);

		currCluster.optimizeReactions();
	}

	return;
}

void NEClusterReactionNetwork::reinitializeConnectivities() {
	// Loop on all the reactants to reset their connectivities
	std::for_each(allReactants.begin(), allReactants.end(),
			[](IReactant& currReactant) {
				currReactant.resetConnectivities();
			});

	return;
}

void NEClusterReactionNetwork::updateConcentrationsFromArray(
		double * concentrations) {

	// Set the concentration on each reactant.
	std::for_each(allReactants.begin(), allReactants.end(),
			[&concentrations](IReactant& currReactant) {
				auto id = currReactant.getId() - 1;
				currReactant.setConcentration(concentrations[id]);
			});

	// Set the moments
	auto const& superTypeMap = getAll(ReactantType::NESuper);
	std::for_each(superTypeMap.begin(), superTypeMap.end(),
			[&concentrations](const ReactantMap::value_type& currMapItem) {
				auto& cluster = static_cast<NESuperCluster&>(*(currMapItem.second));
				cluster.setZerothMoment(concentrations[cluster.getId() - 1]);
				cluster.setMoment(concentrations[cluster.getMomentId() - 1]);
			});

	return;
}

std::vector<std::vector<int> > NEClusterReactionNetwork::getCompositionList() const {
	// Create the list that will be returned
	std::vector<std::vector<int> > compList;

	// Loop on all the reactants
	std::for_each(allReactants.begin(), allReactants.end(),
			[&compList](IReactant& currReactant) {
				// Get the composition
				auto comp = currReactant.getComposition();
				std::vector <int> compVec;
				compVec.push_back(comp[toCompIdx(Species::Xe)]);

				// Save the composition in the list
				compList.push_back(compVec);
			});

	return compList;
}

void NEClusterReactionNetwork::getDiagonalFill(int *diagFill) {

	// Degrees of freedom is the total number of clusters in the network
	const int dof = getDOF();

	// Get the connectivity for each reactant
	std::for_each(allReactants.begin(), allReactants.end(),
			[&diagFill, &dof, this](const IReactant& reactant) {
				// Get the reactant and its connectivity
				auto const& connectivity = reactant.getConnectivity();
				auto connectivityLength = connectivity.size();
				// Get the reactant id so that the connectivity can be lined up in
				// the proper column
				auto id = reactant.getId() - 1;
				// Create the vector that will be inserted into the dFill map
				std::vector<int> columnIds;
				// Add it to the diagonal fill block
				for (int j = 0; j < connectivityLength; j++) {
					// The id starts at j*connectivity length and is always offset
					// by the id, which denotes the exact column.
					auto index = id * dof + j;
					diagFill[index] = connectivity[j];
					// Add a column id if the connectivity is equal to 1.
					if (connectivity[j] == 1) {
						columnIds.push_back(j);
					}
				}
				// Update the map
				dFillMap[id] = columnIds;
			});

	// Get the connectivity for each moment
	for (auto const& currMapItem : getAll(ReactantType::NESuper)) {

		// Get the reactant and its connectivity
		auto const& reactant =
				static_cast<NESuperCluster&>(*(currMapItem.second));

		auto const& connectivity = reactant.getConnectivity();
		auto connectivityLength = connectivity.size();
		// Get the xenon moment id so that the connectivity can be lined up in
		// the proper column
		auto id = reactant.getMomentId() - 1;

		// Create the vector that will be inserted into the dFill map
		std::vector<int> columnIds;
		// Add it to the diagonal fill block
		for (int j = 0; j < connectivityLength; j++) {
			// The id starts at j*connectivity length and is always offset
			// by the id, which denotes the exact column.
			auto index = (id) * dof + j;
			diagFill[index] = connectivity[j];
			// Add a column id if the connectivity is equal to 1.
			if (connectivity[j] == 1) {
				columnIds.push_back(j);
			}
		}
		// Update the map
		dFillMap[id] = columnIds;
	}

	return;
}

void NEClusterReactionNetwork::computeRateConstants() {
	// Local declarations
	double rate = 0.0;
	// Initialize the value for the biggest production rate
	double biggestProductionRate = 0.0;

	// Loop on all the production reactions
	for (auto& currReactionInfo : productionReactionMap) {

		auto& currReaction = currReactionInfo.second;

		// Compute the rate
		rate = calculateReactionRateConstant(*currReaction);
		// Set it in the reaction
		currReaction->kConstant = rate;

		// Check if the rate is the biggest one up to now
		if (rate > biggestProductionRate)
			biggestProductionRate = rate;
	}

	// Loop on all the dissociation reactions
	for (auto& currReactionInfo : dissociationReactionMap) {

		auto& currReaction = currReactionInfo.second;

		// Compute the rate
		rate = calculateDissociationConstant(*currReaction);
		// Set it in the reaction
		currReaction->kConstant = rate;
	}

	// Set the biggest rate
	biggestRate = biggestProductionRate;

	return;
}

void NEClusterReactionNetwork::computeAllFluxes(double *updatedConcOffset) {

	// ----- Compute all of the new fluxes -----
	std::for_each(allReactants.begin(), allReactants.end(),
			[&updatedConcOffset](IReactant& cluster) {

				// Compute the flux
				auto flux = cluster.getTotalFlux();
				// Update the concentration of the cluster
				auto reactantIndex = cluster.getId() - 1;
				updatedConcOffset[reactantIndex] += flux;
			});

	// ---- Moments ----
	for (auto const& currMapItem : getAll(ReactantType::NESuper)) {

		auto& superCluster = static_cast<NESuperCluster&>(*(currMapItem.second));

		// Compute the xenon moment flux
		auto flux = superCluster.getMomentFlux();
		// Update the concentration of the cluster
		auto reactantIndex = superCluster.getMomentId() - 1;
		updatedConcOffset[reactantIndex] += flux;
	}

	return;
}

void NEClusterReactionNetwork::computeAllPartials(double *vals, int *indices,
		int *size) {

	// Initial declarations
	const int dof = getDOF();
	std::vector<double> clusterPartials;
	clusterPartials.resize(dof, 0.0);

	// Update the column in the Jacobian that represents each normal reactant
	for (auto const& superMapItem : getAll(ReactantType::Xe)) {
		auto const& reactant = *(superMapItem.second);

		// Get the reactant index
		auto reactantIndex = reactant.getId() - 1;

		// Get the partial derivatives
		reactant.getPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		auto pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		auto pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];

			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}
	}

	// Get the super clusters
	auto const& superClusters = getAll(ReactantType::NESuper);

	// Update the column in the Jacobian that represents the moment for the super clusters
	for (auto const& currMapItem : superClusters) {
		auto const& reactant =
				static_cast<NESuperCluster&>(*(currMapItem.second));

		// Get the super cluster index
		auto reactantIndex = reactant.getId() - 1;

		// Get the partial derivatives
		reactant.getPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		auto pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		auto pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];
			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}

		// Get the helium moment index
		reactantIndex = reactant.getMomentId() - 1;

		// Get the partial derivatives
		reactant.getMomentPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];
			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}
	}

	return;
}

