#include "NEClusterReactionNetwork.h"
#include "NESuperCluster.h"
#include <xolotlPerf.h>
#include <sstream>
#include <Constants.h>

using namespace xolotlCore;

NEClusterReactionNetwork::NEClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork( { ReactantType::V, ReactantType::I, ReactantType::Xe,
				ReactantType::XeV, ReactantType::XeI, ReactantType::NESuper },
				registry) {

	// Initialize default properties
	dissociationsEnabled = true;

	return;
}

double NEClusterReactionNetwork::calculateDissociationConstant(
		const DissociationReaction& reaction, int i) {

	// If the dissociations are not allowed
	if (!dissociationsEnabled)
		return 0.0;

	// Compute the atomic volume
	double atomicVolume = 0.5 * xolotlCore::uraniumDioxydeLatticeConstant
			* xolotlCore::uraniumDioxydeLatticeConstant
			* xolotlCore::uraniumDioxydeLatticeConstant;

	// Get the rate constant from the reverse reaction
	double kPlus = reaction.reverseReaction->kConstant[i];

	// Calculate and return
	double bindingEnergy = computeBindingEnergy(reaction);
	double k_minus_exp = exp(
			-1.0 * bindingEnergy / (xolotlCore::kBoltzmann * temperature)); // We can use the network temperature
	// because this method is called only
	// when the temperature is updated
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
	auto& singleXeCluster = static_cast<NECluster&>(*(get(Species::Xe, firstSize)));
	// Consider each Xe super cluster.
	for (auto const& currMapItem : getAll(ReactantType::Xe)) {

		auto& xeReactant = static_cast<NECluster&>(*(currMapItem.second));
		// Consider each potential product in normal clusters
		for (auto const& currMapItemBis : getAll(ReactantType::Xe)) {

			auto& product = static_cast<NECluster&>(*(currMapItemBis.second));
			// Check that the reaction can occur
			if (singleXeCluster.getDiffusionFactor() > 0.0
					|| xeReactant.getDiffusionFactor() > 0.0) {
				if (checkOverlap(singleXeCluster, xeReactant, product)) {
					// Create the reaction
					std::unique_ptr<ProductionReaction> reaction(
							new ProductionReaction(singleXeCluster,
									xeReactant));
					auto& prref = add(std::move(reaction));
					// Tell the reactants that they are in this reaction
					singleXeCluster.participateIn(prref, product);
					xeReactant.participateIn(prref, product);
					product.resultFrom(prref, product);

					// Check if the reverse reaction is allowed
					checkForDissociation(&product, prref);
				}
			}
		}
		// Consider each potential product in super clusters
		for (auto const& currMapItemBis : getAll(ReactantType::NESuper)) {

			auto& product = static_cast<NECluster&>(*(currMapItemBis.second));
			// Check that the reaction can occur
			if (singleXeCluster.getDiffusionFactor() > 0.0
					|| xeReactant.getDiffusionFactor() > 0.0) {
				if (checkOverlap(singleXeCluster, xeReactant, product)) {
					// Create the reaction
					std::unique_ptr<ProductionReaction> reaction(
							new ProductionReaction(singleXeCluster,
									xeReactant));
					auto& prref = add(std::move(reaction));
					// Tell the reactants that they are in this reaction
					singleXeCluster.participateIn(prref, product);
					xeReactant.participateIn(prref, product);
					product.resultFrom(prref, product);

					// Check if the reverse reaction is allowed
					checkForDissociation(&product, prref);
				}
			}
		}
	}
	// Consider each Xe super cluster.
	for (auto const& currMapItem : getAll(ReactantType::NESuper)) {

		auto& xeReactant = static_cast<NECluster&>(*(currMapItem.second));
		// Consider each potential product
		for (auto const& currMapItemBis : getAll(ReactantType::NESuper)) {

			auto& product = static_cast<NECluster&>(*(currMapItemBis.second));
			// Check that the reaction can occur
			if (singleXeCluster.getDiffusionFactor() > 0.0
					|| xeReactant.getDiffusionFactor() > 0.0) {
				if (checkOverlap(singleXeCluster, xeReactant, product)) {
					// Create the reaction
					std::unique_ptr<ProductionReaction> reaction(
							new ProductionReaction(singleXeCluster,
									xeReactant));
					auto& prref = add(std::move(reaction));
					// Tell the reactants that they are in this reaction
					singleXeCluster.participateIn(prref, product);
					xeReactant.participateIn(prref, product);
					product.resultFrom(prref, product);

					// Check if the reverse reaction is allowed
					checkForDissociation(&product, prref);
				}
			}
		}
	}

	return;
}

void NEClusterReactionNetwork::checkForDissociation(
		IReactant * emittingReactant, ProductionReaction& reaction) {
	// Check if at least one of the potentially emitted cluster is size one
	if (reaction.first.getSize() != 1 && reaction.second.getSize() != 1) {
		// Don't add the reverse reaction
		return;
	}

	// The reaction can occur, create the dissociation
	// Create a dissociation reaction
	// TODO can this be on the stack?
	std::unique_ptr<DissociationReaction> dissociationReaction(
			new DissociationReaction(*emittingReactant, reaction.first,
					reaction.second));
	// Set the reverse reaction
	dissociationReaction->reverseReaction = &reaction;
	auto& drref = add(std::move(dissociationReaction));
	// Tell the reactants that their are in this reaction
	reaction.first.participateIn(drref, *emittingReactant);
	reaction.second.participateIn(drref, *emittingReactant);
	emittingReactant->emitFrom(drref, *emittingReactant);

	return;
}

void NEClusterReactionNetwork::setTemperature(double temp, int i) {
	ReactionNetwork::setTemperature(temp, i);

	computeRateConstants(i);

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
			});

	// Get all the super clusters and loop on them
	// Have to use allReactants again to be sure the ordering is the same across plateforms
	std::for_each(allReactants.begin(), allReactants.end(),
			[&id, this](IReactant& currReactant) {

				if (currReactant.getType() == ReactantType::NESuper) {
					auto& currCluster = static_cast<NESuperCluster&>(currReactant);
					id++;
					currCluster.setMomentId(id);

					// Update the size
					IReactant::SizeType clusterSize = (double)currCluster.getAverage()
					+ (double)(currCluster.getNTot() - 1) / 2.0;
					if (clusterSize > maxClusterSizeMap[ReactantType::Xe]) {
						maxClusterSizeMap[ReactantType::Xe] = clusterSize;
					}
				}
			});

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

IReactant * NEClusterReactionNetwork::getSuperFromComp(IReactant::SizeType nXe,
		IReactant::SizeType nD, IReactant::SizeType nT,
		IReactant::SizeType nV) const {

	// Requests for finding a particular supercluster have high locality.
	// See if the last supercluster we were asked to find is the right
	// one for this request.
	static IReactant* lastRet;
	if (lastRet and static_cast<NESuperCluster*>(lastRet)->isIn(nXe)) {
		return lastRet;
	}

	// We didn't find the last supercluster in our cache, so do a full lookup.
	IReactant* ret = nullptr;

	for (auto const& superMapItem : getAll(ReactantType::NESuper)) {

		auto const& reactant =
				static_cast<NESuperCluster&>(*(superMapItem.second));
		if (reactant.isIn(nXe)) {
			lastRet = superMapItem.second.get();
			return superMapItem.second.get();
		}
	}

	return ret;
}

void NEClusterReactionNetwork::getDiagonalFill(SparseFillMap& fillMap) {

	// Degrees of freedom is the total number of clusters in the network
	const int dof = getDOF();

	// Get the connectivity for each reactant
	std::for_each(allReactants.begin(), allReactants.end(),
			[&fillMap, &dof, this](const IReactant& reactant) {
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
					// Add a column id if the connectivity is equal to 1.
					if (connectivity[j] == 1) {
						fillMap[id].emplace_back(j);
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
			// Add a column id if the connectivity is equal to 1.
			if (connectivity[j] == 1) {
				fillMap[id].emplace_back(j);
				columnIds.push_back(j);
			}
		}
		// Update the map
		dFillMap[id] = columnIds;
	}

	return;
}

void NEClusterReactionNetwork::computeAllFluxes(double *updatedConcOffset,
		int i) {

	// ----- Compute all of the new fluxes -----
	std::for_each(allReactants.begin(), allReactants.end(),
			[&updatedConcOffset,&i](IReactant& cluster) {

				// Compute the flux
				auto flux = cluster.getTotalFlux(i);
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

void NEClusterReactionNetwork::computeAllPartials(
		const std::vector<size_t>& startingIdx, const std::vector<int>& indices,
		std::vector<double>& vals, int i) const {

	// Initial declarations
	const int dof = getDOF();
	std::vector<double> clusterPartials(dof, 0.0);

	// Update the column in the Jacobian that represents each normal reactant
	for (auto const& superMapItem : getAll(ReactantType::Xe)) {
		auto const& reactant = *(superMapItem.second);

		// Get the reactant index
		auto reactantIndex = reactant.getId() - 1;

		// Get the partial derivatives
		reactant.getPartialDerivatives(clusterPartials, i);
		// Get the list of column ids from the map
		auto const& pdColIdsVector = dFillMap.at(reactantIndex);

		// Loop over the list of column ids
		auto myStartingIdx = startingIdx[reactantIndex];
		for (int j = 0; j < pdColIdsVector.size(); j++) {
			// Get the partial derivative from the array of all of the partials
			vals[myStartingIdx + j] = clusterPartials[pdColIdsVector[j]];

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
		reactant.getPartialDerivatives(clusterPartials, i);

		{
			// Get the list of column ids from the map
			auto const& pdColIdsVector = dFillMap.at(reactantIndex);

			// Loop over the list of column ids
			auto myStartingIdx = startingIdx[reactantIndex];
			for (int j = 0; j < pdColIdsVector.size(); j++) {
				// Get the partial derivative from the array of all of the partials
				vals[myStartingIdx + j] = clusterPartials[pdColIdsVector[j]];

				// Reset the cluster partial value to zero. This is much faster
				// than using memset.
				clusterPartials[pdColIdsVector[j]] = 0.0;
			}
		}
		{
			// Get the Xe momentum index
			auto reactantIndex = reactant.getMomentId() - 1;

			// Get the partial derivatives
			reactant.getMomentPartialDerivatives(clusterPartials);
			// Get the list of column ids from the map
			auto const& pdColIdsVector = dFillMap.at(reactantIndex);

			// Loop over the list of column ids
			auto myStartingIdx = startingIdx[reactantIndex];
			for (int j = 0; j < pdColIdsVector.size(); j++) {
				// Get the partial derivative from the array of all of the partials
				vals[myStartingIdx + j] = clusterPartials[pdColIdsVector[j]];

				// Reset the cluster partial value to zero. This is much faster
				// than using memset.
				clusterPartials[pdColIdsVector[j]] = 0.0;
			}
		}
	}

	return;
}

