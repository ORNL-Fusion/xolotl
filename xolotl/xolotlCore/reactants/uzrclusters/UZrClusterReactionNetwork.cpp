#include <cassert>
#include <iterator>
#include "UZrClusterReactionNetwork.h"
#include "UZrCluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

namespace xolotlCore {

UZrClusterReactionNetwork::UZrClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork(
				{ ReactantType::Xe, ReactantType::V, ReactantType::XeV },
				registry) {

	// Initialize default properties
	dissociationsEnabled = true;

	return;
}

double UZrClusterReactionNetwork::calculateDissociationConstant(
		const DissociationReaction &reaction, int i) {

	// If the dissociations are not allowed
	if (!dissociationsEnabled)
		return 0.0;

	// TODO: Compute the atomic volume correctly depending on the structure of the material
	// You can look at NEClusterReactionNetwork
	double atomicVolume = 0.5 * pow(latticeParameter, 3);

	// Get the rate constant from the reverse reaction
	double kPlus = reaction.reverseReaction->kConstant[i];

	// Calculate and return
	double bindingEnergy = computeBindingEnergy(reaction);
	double k_minus_exp = exp(
			-1.0 * bindingEnergy / (xolotlCore::kBoltzmann * temperature));
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;

	return k_minus;
}

void UZrClusterReactionNetwork::createReactionConnectivity() {
	// TODO: all the reactions are defined here
	// I put simple clustering as an example

	// Initial declarations
	IReactant::SizeType firstSize = 0, secondSize = 0, productSize = 0;

	// Single species clustering (Xe, V)
	// X_a + X_i --> X_(a+i)
	// Make a vector of types
	std::vector<ReactantType> typeVec { ReactantType::Xe, ReactantType::V };
	// Loop on it
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {

		auto currType = *tvIter;

		// Consider all reactants of the current type.
		auto const &currTypeReactantMap = getAll(currType);
		for (auto firstIt = currTypeReactantMap.begin();
				firstIt != currTypeReactantMap.end(); firstIt++) {

			auto &firstReactant = *(firstIt->second);
			// Get its size
			firstSize = firstReactant.getSize();

			// Loop on the second cluster starting at the same pointer to avoid double counting
			for (auto secondIt = firstIt; secondIt != currTypeReactantMap.end();
					secondIt++) {

				auto &secondReactant = *(secondIt->second);

				// At least one of them should diffuse
				if (firstReactant.getDiffusionFactor() <= 0.0
						&& secondReactant.getDiffusionFactor() <= 0.0)
					continue;

				// Get its size
				secondSize = secondReactant.getSize();
				productSize = firstSize + secondSize;

				// Get the product
				auto product = get(toSpecies(currType), productSize);
				// Check that the reaction can occur
				if (product) {
					// Also checks on dissociation
					defineProductionReaction(firstReactant, secondReactant,
							*product);
				}
			}
		}
	}

	// Xenon-Vacancy clustering
	// Xe_a + V_b --> (Xe_a)(V_b)
	// Consider each Xe reactant.
	for (auto const &xeMapItem : getAll(ReactantType::Xe)) {

		auto &xeReactant = *(xeMapItem.second);
		// Get its size
		firstSize = xeReactant.getSize();

		// Consider product with each V cluster
		for (auto const &vMapItem : getAll(ReactantType::V)) {

			auto &vReactant = *(vMapItem.second);

			// At least one of them should diffuse
			if (xeReactant.getDiffusionFactor() <= 0.0
					&& vReactant.getDiffusionFactor() <= 0.0)
				continue;

			// Get its size
			secondSize = vReactant.getSize();

			// Check if product exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::Xe)] = firstSize;
			newComp[toCompIdx(Species::V)] = secondSize;
			auto product = get(ReactantType::XeV, newComp);

			// Create the reaction
			if (product) {
				// Also checks on dissociation
				defineProductionReaction(xeReactant, vReactant, *product);
			}
		}
	}

	// Vacancy absorption by XeV clusters
	// (Xe_a)(V_b) + V_c --> (Xe_a)[V_(b+c)]
	// Consider each V cluster.
	for (auto const &vMapItem : getAll(ReactantType::V)) {

		auto &vReactant = *(vMapItem.second);

		// Skip if it can't diffuse because mixed clusters can't diffuse
		if (xolotlCore::equal(vReactant.getDiffusionFactor(), 0.0))
			continue;
		// Get the V size
		firstSize = vReactant.getSize();
		// Consider product with every XeV cluster.
		for (auto const &xeVMapItem : getAll(ReactantType::XeV)) {

			auto &xeVReactant = *(xeVMapItem.second);

			// Get its composition
			auto &comp = xeVReactant.getComposition();
			// Create the composition of the potential product
			int newNumXe = comp[toCompIdx(Species::Xe)];
			int newNumV = comp[toCompIdx(Species::V)] + firstSize;

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::Xe)] = newNumXe;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::XeV, newComp);

			// Create the reaction
			if (product) {
				// Also checks on dissociation
				defineProductionReaction(vReactant, xeVReactant, *product);
			}
		}
	}

	// Xenon absorption by XeV clusters
	// (Xe_a)(V_b) + Xe_c --> [Xe_(a+c)](V_b)
	// Consider each Xe cluster.
	for (auto const &xeMapItem : getAll(ReactantType::Xe)) {

		auto &xeReactant = *(xeMapItem.second);

		// Skip if it can't diffuse because mixed clusters can't diffuse
		if (xolotlCore::equal(xeReactant.getDiffusionFactor(), 0.0))
			continue;
		// Get the Xe size
		firstSize = xeReactant.getSize();
		// Consider product with every XeV cluster.
		for (auto const &xeVMapItem : getAll(ReactantType::XeV)) {

			auto &xeVReactant = *(xeVMapItem.second);

			// Get its composition
			auto &comp = xeVReactant.getComposition();
			// Create the composition of the potential product
			int newNumXe = comp[toCompIdx(Species::Xe)] + firstSize;
			int newNumV = comp[toCompIdx(Species::V)];

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::Xe)] = newNumXe;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::XeV, newComp);

			// Create the reaction
			if (product) {
				// Also checks on dissociation
				defineProductionReaction(xeReactant, xeVReactant, *product);
			}
		}
	}

	return;
}

void UZrClusterReactionNetwork::setTemperature(double temp, int i) {
	ReactionNetwork::setTemperature(temp, i);

	computeRateConstants(i);

	return;
}

void UZrClusterReactionNetwork::reinitializeNetwork() {

	// Reset the Ids
	// std::for_each is guaranteed to visit reactants in order for C++11.
	int id = 0;
	std::for_each(allReactants.begin(), allReactants.end(),
			[&id](IReactant &currReactant) {
				id++;
				currReactant.setId(id);
			});

	return;
}

void UZrClusterReactionNetwork::reinitializeConnectivities() {

	// Reset connectivities of each reactant.
	std::for_each(allReactants.begin(), allReactants.end(),
			[](IReactant &currReactant) {
				currReactant.resetConnectivities();
			});

	return;
}

void UZrClusterReactionNetwork::updateConcentrationsFromArray(
		double *concentrations) {

	// Set the concentration on each reactant.
	std::for_each(allReactants.begin(), allReactants.end(),
			[&concentrations](IReactant &currReactant) {
				auto id = currReactant.getId() - 1;
				currReactant.setConcentration(concentrations[id]);
			});

	return;
}

std::vector<std::vector<int> > UZrClusterReactionNetwork::getCompositionList() const {
	// Create the list that will be returned
	std::vector<std::vector<int> > compList;

	// Loop on all the reactants
	std::for_each(allReactants.begin(), allReactants.end(),
			[&compList](IReactant &currReactant) {
				// Get the composition
				auto comp = currReactant.getComposition();
				std::vector<int> compVec;
				compVec.push_back(comp[toCompIdx(Species::Xe)]);
				compVec.push_back(comp[toCompIdx(Species::V)]);

				// Save the composition in the list
				compList.push_back(compVec);
			});

	return compList;
}

void UZrClusterReactionNetwork::getDiagonalFill(SparseFillMap &fillMap) {
	// Degrees of freedom is the total number of clusters in the network
	const int dof = getDOF();

	// Get the connectivity for each reactant
	std::for_each(allReactants.begin(), allReactants.end(),
			[&fillMap, &dof, this](const IReactant &reactant) {

				// Get the reactant's connectivity
				auto const &connectivity = reactant.getConnectivity();
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

	return;
}

double UZrClusterReactionNetwork::getTotalAtomConcentration(int i) {
	// Initial declarations
	double conc = 0.0;

	// Sum over all Xe clusters.
	for (auto const &currMapItem : getAll(ReactantType::Xe)) {

		// Get the cluster and its composition
		auto const &cluster = *(currMapItem.second);
		double size = cluster.getSize();

		// Add the concentration times the Xe content to the total xenon concentration
		conc += cluster.getConcentration() * size;
	}

	// Sum over all XeV clusters.
	for (auto const &currMapItem : getAll(ReactantType::XeV)) {

		// Get the cluster and its composition
		auto const &cluster = *(currMapItem.second);
		auto &comp = cluster.getComposition();

		// Add the concentration times the Xe content to the total xenon concentration
		conc += cluster.getConcentration() * comp[toCompIdx(Species::Xe)];
	}

	return conc;
}

double UZrClusterReactionNetwork::getTotalTrappedAtomConcentration(int i) {
	// Initial declarations
	double conc = 0.0;

	// Sum over all XeV clusters.
	for (auto const &currMapItem : getAll(ReactantType::XeV)) {
		// Get the cluster and its composition
		auto const &cluster = *(currMapItem.second);
		auto &comp = cluster.getComposition();

		// Add the concentration times the Xe content to the total xenon concentration
		conc += cluster.getConcentration() * comp[toCompIdx(Species::Xe)];
	}

	return conc;
}

double UZrClusterReactionNetwork::getTotalVConcentration() {
	// Initial declarations
	double vConc = 0.0;

	// Sum over all V clusters.
	for (auto const &currMapItem : getAll(ReactantType::V)) {
		// Get the cluster and its composition
		auto const &cluster = *(currMapItem.second);
		double size = cluster.getSize();

		// Add the concentration times the V content to the total vacancy concentration
		vConc += cluster.getConcentration() * size;
	}

	// Sum over all XeV clusters
	for (auto const &currMapItem : getAll(ReactantType::XeV)) {
		// Get the cluster and its composition
		auto const &cluster = *(currMapItem.second);
		auto &comp = cluster.getComposition();

		// Add the concentration times the V content to the total vacancy concentration
		vConc += cluster.getConcentration() * comp[toCompIdx(Species::V)];
	}

	return vConc;
}

double UZrClusterReactionNetwork::getTotalIConcentration() {
	return 0.0;
}

void UZrClusterReactionNetwork::computeAllFluxes(double *updatedConcOffset,
		int i) {

	// ----- Compute all of the new fluxes -----
	std::for_each(allReactants.begin(), allReactants.end(),
			[&updatedConcOffset, &i](IReactant &cluster) {
				// Compute the flux
				auto flux = cluster.getTotalFlux(i);
				// Update the concentration of the cluster
				auto reactantIndex = cluster.getId() - 1;
				updatedConcOffset[reactantIndex] += flux;
			});

	return;
}

void UZrClusterReactionNetwork::computeAllPartials(
		const std::vector<size_t> &startingIdx, const std::vector<int> &indices,
		std::vector<double> &vals, int i) const {
	// Initial declarations
	const int dof = getDOF();
	std::vector<double> clusterPartials(dof, 0.0);

	// Make a vector of types for the non super clusters
	std::vector<ReactantType> typeVec { ReactantType::Xe, ReactantType::V,
			ReactantType::XeV };
	// Loop on it
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {

		auto currType = *tvIter;

		// Consider all reactants of the current type.
		auto const &currTypeReactantMap = getAll(currType);

		// Update the column in the Jacobian that represents each normal reactant
		for (auto const &currMapItem : currTypeReactantMap) {

			auto const &reactant =
					static_cast<UZrCluster&>(*(currMapItem.second));

			// Get the reactant index
			auto reactantIndex = reactant.getId() - 1;

			// Get the partial derivatives
			reactant.getPartialDerivatives(clusterPartials, i);
			// Get the list of column ids from the map
			auto const &pdColIdsVector = dFillMap.at(reactantIndex);

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

double UZrClusterReactionNetwork::computeBindingEnergy(
		const DissociationReaction &reaction) const {

	// TODO: change to the desired formula or base it on formation energies
	double bindingEnergy = 5.0;
	if (reaction.dissociating.getType() == ReactantType::Xe
			&& reaction.first.getType() == ReactantType::Xe) {
		if (reaction.dissociating.getSize() == 2)
			bindingEnergy = 0.5;
		else
			bindingEnergy = 1.0;
	}
	if (reaction.dissociating.getType() == ReactantType::V
			&& reaction.first.getType() == ReactantType::V) {
		int size = reaction.dissociating.getSize();
		bindingEnergy = 1.73
				- 2.59
						* (pow((double) size, 2.0 / 3.0)
								- pow((double) size - 1.0, 2.0 / 3.0));
	}
	if ((reaction.dissociating.getType() == ReactantType::XeV)
			&& (reaction.first.getType() == ReactantType::V
					|| reaction.second.getType() == ReactantType::V)) {
		auto &comp = reaction.dissociating.getComposition();
		bindingEnergy = 1.73
				- 2.59
						* (pow((double) comp[toCompIdx(Species::V)], 2.0 / 3.0)
								- pow(
										(double) comp[toCompIdx(Species::V)]
												- 1.0, 2.0 / 3.0))
				+ 2.5
						* log(
								1.0
										+ ((double) comp[toCompIdx(Species::Xe)]
												/ (double) comp[toCompIdx(
														Species::V)]));
	}

//	if (bindingEnergy < -5.0)
//	std::cout << "dissociation: " << reaction.dissociating.getName() << " -> "
//			<< reaction.first.getName() << " + "
//			<< reaction.second.getName() << " : " << bindingEnergy
//			<< std::endl;

	return max(bindingEnergy, -5.0);
}

} // namespace xolotlCore

