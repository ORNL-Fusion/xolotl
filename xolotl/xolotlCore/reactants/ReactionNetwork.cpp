#include "ReactionNetwork.h"
#include <xolotlPerf.h>
#include <iostream>
#include <cassert>

namespace xolotlCore {

ReactionNetwork::ReactionNetwork(
		const std::set<ReactantType>& _knownReactantTypes,
		ReactantType _superClusterType,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> _registry) :
		knownReactantTypes(_knownReactantTypes), superClusterType(
				_superClusterType), handlerRegistry(_registry), temperature(
				0.0), dissociationsEnabled(true) {

	// Ensure our per-type cluster map can store Reactants of the types
	// we support.
	for (auto const& currType : knownReactantTypes) {
		clusterTypeMap.emplace(
				std::make_pair(currType, IReactionNetwork::ReactantMap()));
	}

	// Ensure we have a baseline for determining max cluster size for
	// the types we support.
	for (auto const& currType : knownReactantTypes) {
		maxClusterSizeMap.insert( { currType, 0 });
	}
	return;
}

void ReactionNetwork::add(std::unique_ptr<IReactant> reactant) {

	// Ensure we have a valid object to work with.
	assert(reactant);

	// Get the composition
	auto& composition = reactant->getComposition();

	// Check if we already know about this reactant.
	auto& currTypeMap = clusterTypeMap[reactant->getType()];
	auto iter = currTypeMap.find(composition);
	if (iter == currTypeMap.end()) {

		// Set the id for this cluster
		// (It is networkSize+1 because we haven't added
		// it to the network yet.)
		reactant->setId(size() + 1);

		// Update the max cluster size for the new cluster's type.
		maxClusterSizeMap[reactant->getType()] = std::max(reactant->getSize(),
				maxClusterSizeMap[reactant->getType()]);

		// Note the reactant in our flat list of all reactants.
		allReactants.emplace_back(*reactant);

		// Give reactant to the appropriate per-type map.
		currTypeMap.emplace(composition, std::move(reactant));

	} else {
		std::stringstream errStream;
		errStream << "ReactionNetwork: not adding duplicate "
				<< toString(reactant->getType()) << ':' << composition;

		throw errStream.str();
	}

	return;
}

// TODO have this return an iterator to the appropriate map.  Probably
// will have to have it return a pair with iter and bool, to ease
// detection of whether we found it or not.  (Figuring out which
// map to test the end against is probable less easy.)
// Alternative is to remove 'get()' entirely and rely on callers
// to getAll() for particular type, then do their find/test to see if
// we found it.  That would be (slightly?) faster since now caller
// has to test return for nullptr anyway, and we're already testing
// against end() of the appropriate map.
IReactant * ReactionNetwork::get(ReactantType type,
		const IReactant::Composition& comp) const {

	IReactant* ret = nullptr;

	// Check if the reactant is in the map
	auto const& currTypeMap = clusterTypeMap.at(type);
	auto iter = currTypeMap.find(comp);
	if (iter != currTypeMap.end()) {
		ret = iter->second.get();
	}

	return ret;
}

double ReactionNetwork::calculateReactionRateConstant(
		const ProductionReaction& reaction, int i) const {
	// Get the reaction radii
	double r_first = reaction.first.getReactionRadius();
	double r_second = reaction.second.getReactionRadius();

	// Get the diffusion coefficients
	double firstDiffusion = reaction.first.getDiffusionCoefficient(i);
	double secondDiffusion = reaction.second.getDiffusionCoefficient(i);

	// Calculate and return
	double k_plus = 4.0 * xolotlCore::pi * (r_first + r_second)
			* (firstDiffusion + secondDiffusion);

	return k_plus;
}

void ReactionNetwork::fillConcentrationsArray(double * concentrations) {

	// Fill the array
	std::for_each(allReactants.begin(), allReactants.end(),
			[&concentrations](const IReactant& currReactant) {
				auto id = currReactant.getId() - 1;
				concentrations[id] = currReactant.getConcentration();
			});

	return;
}

void ReactionNetwork::updateConcentrationsFromArray(double * concentrations) {

	std::for_each(allReactants.begin(), allReactants.end(),
			[&concentrations](IReactant& currReactant) {
				auto id = currReactant.getId() - 1;
				currReactant.setConcentration(concentrations[id]);
			});

	return;
}

void ReactionNetwork::setTemperature(double temp, int i) {
	// Set the temperature
	temperature = temp;

	// Update the temperature for all of the clusters
	std::for_each(allReactants.begin(), allReactants.end(),
			[&temp,&i](IReactant& currReactant) {

				// This part will set the temperature in each reactant
				// and recompute the diffusion coefficient
				currReactant.setTemperature(temp, i);
			});

	return;
}

ProductionReaction& ReactionNetwork::add(
		std::unique_ptr<ProductionReaction> reaction) {
	// Ensure we know about the reaction.
	// Map's emplace() returns a pair (iter, bool) where
	// iter points to the item in the map and the bool indicates
	// whether it was added by this emplace() call.
	auto key = reaction->descriptiveKey();
	auto eret = productionReactionMap.emplace(key, std::move(reaction));
	// Regardless of whether we added it in this emplace() call or not,
	// the iter within eret refers to the desired reaction in the map.
	return *(eret.first->second);
}

DissociationReaction& ReactionNetwork::add(
		std::unique_ptr<DissociationReaction> reaction) {
	// Unlike addProductionReaction, we use a check-add approach
	// instead of emplace() because if the item isn't already in
	// the network, we will need to contruct its reverse reaction.

	// Check if we already know about this reaction.
	auto key = reaction->descriptiveKey();
	auto iter = dissociationReactionMap.find(key);
	if (iter != dissociationReactionMap.end()) {
		// We already knew about the reaction.
		// Return the existing one.
		return *(iter->second);
	}

	// Add the dissociation reaction to our set of known reactions.
	auto eret = dissociationReactionMap.emplace(key, std::move(reaction));

	// Since we checked earlier and the reaction wasn't in the map,
	// our emplace() call should have added it.
	assert(eret.second);

	// Return the newly-added dissociation reaction.
	return *(eret.first->second);
}

void ReactionNetwork::removeReactants(
		const IReactionNetwork::ReactantVector& doomedReactants) {

	// Build a ReactantMatcher functor for the doomed reactants.
	// Doing this here allows us to construct the canonical composition
	// strings for the doomed reactants once and reuse them.
	// If we used an anonymous functor object in the std::remove_if
	// calls we would build these strings several times in this function.
	ReactionNetwork::ReactantMatcher doomedReactantMatcher(doomedReactants);

	// Remove the doomed reactants from the type-specific cluster vectors.
	// First, determine all cluster types used by clusters in the collection
	// of doomed reactants...
	std::set<ReactantType> typesUsed;
	for (IReactant const& reactant : doomedReactants) {
		typesUsed.insert(reactant.getType());
	}

	auto first = allReactants.begin();
	auto last = allReactants.end();
	auto result = first;
	while (first != last) {
		IReactant& reactant = (*first);
		if (!(doomedReactantMatcher(reactant)
				&& typesUsed.find(reactant.getType()) != typesUsed.end())) {
			*result = move(*first);
			++result;
		}
		++first;
	}

	allReactants.erase(result, allReactants.end());

	// ...Next, examine each type's collection of clusters and remove the
	// doomed reactants.
	for (auto currType : typesUsed) {
		auto& clusters = clusterTypeMap[currType];

		for (IReactant const& currDoomedReactant : doomedReactants) {
			auto iter = clusters.find(currDoomedReactant.getComposition());
			assert(iter != clusters.end());
			clusters.erase(iter);
		}
	}

	return;
}

void ReactionNetwork::computeRateConstants(int i) {
	// Local declarations
	double rate = 0.0;
	// Initialize the value for the biggest production rate
	double biggestProductionRate = 0.0;

	// Loop on all the production reactions
	for (auto& currReactionInfo : productionReactionMap) {

		auto& currReaction = currReactionInfo.second;

		// Compute the rate
		rate = calculateReactionRateConstant(*currReaction, i);
		// Set it in the reaction
		currReaction->kConstant[i] = rate;

		// Check if the rate is the biggest one up to now
		if (rate > biggestProductionRate)
			biggestProductionRate = rate;
	}

	// Loop on all the dissociation reactions
	for (auto& currReactionInfo : dissociationReactionMap) {

		auto& currReaction = currReactionInfo.second;

		// Compute the rate
		rate = calculateDissociationConstant(*currReaction, i);

		// Set it in the reaction
		currReaction->kConstant[i] = rate;
	}

	// Set the biggest rate
	biggestRate = biggestProductionRate;

	return;
}

void ReactionNetwork::addGridPoints(int i) {
	// Add grid points to the diffusing clusters first
	for (IReactant& currReactant : allReactants) {
		currReactant.addGridPoints(i);
	}

	// Add grid points
	if (i > 0) {
		while (i > 0) {
			// Loop on all the production reactions
			for (auto& currReactionInfo : productionReactionMap) {
				currReactionInfo.second->kConstant.emplace(
						currReactionInfo.second->kConstant.begin(), 0.0);

			}

			// Loop on all the dissociation reactions
			for (auto& currReactionInfo : dissociationReactionMap) {
				currReactionInfo.second->kConstant.emplace(
						currReactionInfo.second->kConstant.begin(), 0.0);

			}

			// Decrease i
			i--;
		}
	} else {
		// Loop on all the production reactions
		for (auto& currReactionInfo : productionReactionMap) {
			currReactionInfo.second->kConstant.erase(
					currReactionInfo.second->kConstant.begin(),
					currReactionInfo.second->kConstant.begin() - i);
		}
		// Loop on all the dissociation reactions
		for (auto& currReactionInfo : dissociationReactionMap) {
			currReactionInfo.second->kConstant.erase(
					currReactionInfo.second->kConstant.begin(),
					currReactionInfo.second->kConstant.begin() - i);
		}
	}

	return;
}

size_t ReactionNetwork::initPartialsSizes(std::vector<int>& size,
		std::vector<size_t>& startingIdx) const {

	size_t currStartingIdx = 0;

	// Determine the number of items owned by each reactant.
	for (auto idx = 0; idx < getDOF(); ++idx) {
		// Note the starting index of items owned by this reactant.
		startingIdx[idx] = currStartingIdx;

		// Check whether this reactant has any partials.
		auto iter = dFillMap.find(idx);
		if (iter != dFillMap.end()) {
			// The current reactant has some partials values to compute.
			auto currSize = iter->second.size();

			// Note number of valid partial derivatives for this reactant.
			size[idx] = currSize;

			// Advance the starting index past items owned by this reactant.
			currStartingIdx += currSize;
		} else {
			// The current reactant has no partials.
			size[idx] = 0;
		}
	}

	return currStartingIdx;
}

void ReactionNetwork::initPartialsIndices(const std::vector<int>& size,
		const std::vector<size_t>& startingIdx,
		std::vector<int>& indices) const {
	// Determine the number of items owned by each reactant.
	for (auto idx = 0; idx < getDOF(); ++idx) {
		if (size[idx] > 0) {
			// Determine column ids of needed partial derivatives.
			auto const& pdColIdsVector = dFillMap.at(idx);

			// Place the column ids in our location(s) in the indices array.
			auto myStartingIdx = startingIdx[idx];
			for (auto j = 0; j < size[idx]; ++j) {
				indices[myStartingIdx + j] = pdColIdsVector[j];
			}
		}
	}
}

void ReactionNetwork::dumpTo(std::ostream& os) const {
	// Dump flat view of reactants.
	os << size() << " reactants:\n";
	for (auto const& currReactant : allReactants) {
		os << currReactant << '\n';
	}

	// Dump ProductionReactions.
	os << productionReactionMap.size() << " production reactions:\n";
	for (auto const& currMapItem : productionReactionMap) {
		os << *(currMapItem.second) << '\n';
	}

	// Dump DissociationReactions.
	os << dissociationReactionMap.size() << " dissociation reactions:\n";
	for (auto const& currMapItem : dissociationReactionMap) {
		os << *(currMapItem.second) << '\n';
	}

	// For each reactant, dump coefficients it uses for reactions it
	// participates in.
	os << size() << " reactant coefficients:\n";
	for (IReactant const& currReactant : allReactants) {
		currReactant.outputCoefficientsTo(os);
	}
}

} // xolotlCore
