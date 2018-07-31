#include <cassert>
#include <iterator>
#include "PSIClusterReactionNetwork.h"
#include "PSICluster.h"
#include "PSISuperCluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

namespace xolotlCore {

PSIClusterReactionNetwork::PSIClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork( { ReactantType::V, ReactantType::I, ReactantType::He,
				ReactantType::D, ReactantType::T, ReactantType::HeI,
				ReactantType::PSIMixed, ReactantType::PSISuper },
				ReactantType::PSISuper, registry) {

	// Initialize default properties
	dissociationsEnabled = true;

	return;
}

double PSIClusterReactionNetwork::calculateDissociationConstant(
		const DissociationReaction& reaction) const {

	// If the dissociations are not allowed
	if (!dissociationsEnabled)
		return 0.0;

	// The atomic volume is computed by considering the BCC structure of the
	// tungsten. In a given lattice cell in tungsten there are tungsten atoms
	// at each corner and a tungsten atom in the center. The tungsten atoms at
	// the corners are shared across a total of eight cells. The fraction of
	// the volume of the lattice cell that is filled with tungsten atoms is the
	// atomic volume and is a_0^3/(8*1/8 + 1) = 0.5*a_0^3.
	double atomicVolume = 0.5 * xolotlCore::tungstenLatticeConstant
			* xolotlCore::tungstenLatticeConstant
			* xolotlCore::tungstenLatticeConstant;

	// Get the rate constant from the reverse reaction
	double kPlus = reaction.reverseReaction->kConstant;

	// Calculate and return
	double bindingEnergy = computeBindingEnergy(reaction);
	double k_minus_exp = exp(
			-1.0 * bindingEnergy / (xolotlCore::kBoltzmann * temperature));
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;

	return k_minus;
}

void PSIClusterReactionNetwork::defineProductionReactions(IReactant& r1,
		IReactant& r2,
		const std::vector<PendingProductionReactionInfo>& pendingPRInfos,
		bool secondProduct) {

	// Define the production reaction to the network.
	// Do this *once* for the given reactants, since it doesn't
	// depend on the product or the other parameters.
	std::unique_ptr<ProductionReaction> reaction(
			new ProductionReaction(r1, r2));
	auto& prref = add(std::move(reaction));

	// Determine if reverse reaction is allowed.
	auto dissociationAllowed = canDissociate(pendingPRInfos[0].product, prref);

	// Build the product-to-production map that we will
	// use for batched resultsFrom() and defineDissociationReactions() calls.
	ProductToProductionMap prodMap;
	std::for_each(pendingPRInfos.begin(), pendingPRInfos.end(),
			[&prodMap](const PendingProductionReactionInfo& currPRI) {

				prodMap[&(currPRI.product)].emplace_back(currPRI);
			});

	// Tell both reactants they are involved in this reaction
	// with the given parameters if this is not the second product
	if (!secondProduct) {
		prref.first.participateIn(prref, pendingPRInfos);
		prref.second.participateIn(prref, pendingPRInfos);
	}

	// Tell all products they are involved in this reaction
	// with their given parameters.
	std::for_each(prodMap.begin(), prodMap.end(),
			[&prref](const ProductToProductionMap::value_type& prodMapItem) {

				IReactant& currProduct = *(prodMapItem.first);
				ProductToProductionMap::mapped_type const& currPRIs = prodMapItem.second;

				currProduct.resultFrom(prref, currPRIs);
			});

	// Determine if reverse reaction is allowed. Never for the second product
	if (dissociationAllowed && !secondProduct) {
		// Forward reaction can dissociate.
		// Define all dissociation reactions for this forward reaction
		defineDissociationReactions(prref, prodMap);
	}
}

void PSIClusterReactionNetwork::defineAnaProductionReactions(IReactant& r1,
		IReactant& r2, IReactant& product, bool secondProduct) {
	// Define the production reaction to the network.
	// Do this *once* for the given reactants, since it doesn't
	// depend on the product or the other parameters.
	std::unique_ptr<ProductionReaction> reaction(
			new ProductionReaction(r1, r2));
	auto& prref = add(std::move(reaction));

	// Determine if reverse reaction is allowed.
	auto dissociationAllowed = canDissociate(product, prref);

	// Tell both reactants they are involved in this reaction
	// with the given parameters if this is not the second product
	if (!secondProduct) {
		prref.first.participateIn(prref, product);
		prref.second.participateIn(prref, product);
		product.resultFrom(prref, product);
	}

	// Determine if reverse reaction is allowed. Never for the second product
	if (dissociationAllowed && !secondProduct) {
		// Forward reaction can dissociate.
		// Define all dissociation reactions for this forward reaction
		defineAnaDissociationReactions(prref, product);
	}
}

void PSIClusterReactionNetwork::defineDissociationReactions(
		ProductionReaction& forwardReaction,
		const ProductToProductionMap& prodMap) {

	// Consider each production of the given forward reaction.
	std::for_each(prodMap.begin(), prodMap.end(),
			[this,&forwardReaction](const ProductToProductionMap::value_type& currMapItem) {
				// Add a dissociation reaction to our network.
				// Do this once here for each forward reaction product.
				IReactant& emitting = *(currMapItem.first);

				std::unique_ptr<DissociationReaction> dissociationReaction(new DissociationReaction(emitting,
								forwardReaction.first, forwardReaction.second, &forwardReaction));
				auto& drref = add(std::move(dissociationReaction));

				// Tell all participants in this reaction of their involvement.
				ProductToProductionMap::mapped_type const& currPRIs = currMapItem.second;
				drref.first.participateIn(drref, currPRIs);
				drref.second.participateIn(drref, currPRIs);
				emitting.emitFrom(drref, currPRIs);
			});
}

void PSIClusterReactionNetwork::defineAnaDissociationReactions(
		ProductionReaction& forwardReaction, IReactant& emitting) {

	std::unique_ptr<DissociationReaction> dissociationReaction(
			new DissociationReaction(emitting, forwardReaction.first,
					forwardReaction.second, &forwardReaction));
	auto& drref = add(std::move(dissociationReaction));

	// Tell the reactants that they are in this reaction
	forwardReaction.first.participateIn(drref, emitting);
	forwardReaction.second.participateIn(drref, emitting);
	emitting.emitFrom(drref, emitting);
}

void PSIClusterReactionNetwork::createReactionConnectivity() {
	// Initial declarations
	IReactant::SizeType firstSize = 0, secondSize = 0, productSize = 0, maxI =
			getAll(ReactantType::I).size();

	// Single species clustering (He, D, T, V, I)
	// X_(a-i) + X_i --> X_a
	// Make a vector of types
	std::vector<ReactantType> typeVec { ReactantType::He, ReactantType::V,
			ReactantType::I, ReactantType::D, ReactantType::T };
	// Loop on it
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {

		auto currType = *tvIter;

		// Consider all reactants of the current type.
		auto const& currTypeReactantMap = getAll(currType);
		for (auto firstIt = currTypeReactantMap.begin();
				firstIt != currTypeReactantMap.end(); firstIt++) {

			auto& firstReactant = *(firstIt->second);

			// Get its size
			firstSize = firstReactant.getSize();
			// Loop on the second cluster starting at the same pointer to avoid double counting
			for (auto secondIt = firstIt; secondIt != currTypeReactantMap.end();
					secondIt++) {

				auto& secondReactant = *(secondIt->second);

				// Get its size
				secondSize = secondReactant.getSize();
				productSize = firstSize + secondSize;
				// Get the product
				auto product = get(toSpecies(currType), productSize);
				// Check that the reaction can occur
				if (product
						&& (firstReactant.getDiffusionFactor() > 0.0
								|| secondReactant.getDiffusionFactor() > 0.0)) {

					defineProductionReaction(firstReactant, secondReactant,
							*product);
				}
			}
		}
	}

	// Helium absorption by Mixed clusters
	// He_(a) + (He_b)()(V_c) --> [He_(a+b)]()(V_c)
	// Consider each He reactant.
	for (auto const& heMapItem : getAll(ReactantType::He)) {

		auto& heReactant = *(heMapItem.second);

		// Skip if it can't diffuse
		if (xolotlCore::equal(heReactant.getDiffusionFactor(), 0.0))
			continue;
		// Get its size
		firstSize = heReactant.getSize();

		// Consider product with each HeV cluster
		for (auto const& heVMapItem : getAll(ReactantType::PSIMixed)) {

			auto& heVReactant = *(heVMapItem.second);

			// Get its composition
			auto& comp = heVReactant.getComposition();
			// Create the composition of the potential product
			int newNumHe = comp[toCompIdx(Species::He)] + firstSize;
			int newNumD = comp[toCompIdx(Species::D)];
			int newNumT = comp[toCompIdx(Species::T)];
			int newNumV = comp[toCompIdx(Species::V)];

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::D)] = newNumD;
			newComp[toCompIdx(Species::T)] = newNumT;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::PSIMixed, newComp);

			// Check if the product can be a super cluster
			if (!product) {
				// Check if it is a super cluster from the map
				product = getSuperFromComp(newNumHe, newNumD, newNumT, newNumV);
			}
			// Check that the reaction can occur
			if (product
					&& (heReactant.getDiffusionFactor() > 0.0
							|| heVReactant.getDiffusionFactor() > 0.0)) {

				int a[4] = { newNumHe, newNumD, newNumT, newNumV };
				defineProductionReaction(heReactant, heVReactant, *product, a);
			}
		}

		// Consider product with each super cluster
		for (auto const& superMapItem : getAll(ReactantType::PSISuper)) {
			auto& superCluster =
					static_cast<PSISuperCluster&>(*(superMapItem.second));

			// Loop on the potential products
			for (auto const& superMapItemProd : getAll(ReactantType::PSISuper)) {
				auto& superProd =
						static_cast<PSISuperCluster&>(*(superMapItemProd.second));

				// Check if the reactions overlap
				auto& heCluster = static_cast<PSICluster&>(heReactant);
				bool overlap = checkOverlap(heCluster, superCluster, superProd);

				// Skip if they don't overlap
				if (!overlap)
					continue;

				// Check if the super clusters are full
				if (superCluster.isFull() && superProd.isFull()) {
					// This method will check if the reaction is possible and then add it to the list
					defineAnaProductionReactions(heReactant, superCluster,
							superProd);
				} else {
					std::vector<PendingProductionReactionInfo> prInfos;

					// Get the coordinates of the reactant and loop on them
					auto coords = superCluster.getCoordList();
					for (auto const& pair : coords) {
						// Assume the product can only be a super cluster here
						int newNumHe = std::get<0>(pair) + firstSize;
						int newNumD = std::get<1>(pair);
						int newNumT = std::get<2>(pair);
						int newNumV = std::get<3>(pair);
						if (superProd.isIn(newNumHe, newNumD, newNumT, newNumV)
								&& (heReactant.getDiffusionFactor() > 0.0
										|| superCluster.getDiffusionFactor()
												> 0.0)) {
							// Note that current reactant reacts with
							// current superCluster to produce product,
							// according to current parameters.
							int a[4] = { newNumHe, newNumD, newNumT, newNumV };
							int b[4] = { std::get<0>(pair), std::get<1>(pair),
									std::get<2>(pair), std::get<3>(pair) };
							prInfos.emplace_back(superProd, a, b);
						}
					}
					// Now that we know how current reactant reacts with
					// current superCluster, create the production
					// reaction(s) for them.
					if (prInfos.size() > 0)
						defineProductionReactions(heReactant, superCluster,
								prInfos);
				}
			}
		}
	}

	// Vacancy absorption by Mixed clusters
	// (He_a)()(V_b) + V_c --> (He_a)()[V_(b+c)]
	// Consider each V cluster.
	for (auto const& vMapItem : getAll(ReactantType::V)) {

		auto& vReactant = *(vMapItem.second);

		// Skip if it can't diffuse
		if (xolotlCore::equal(vReactant.getDiffusionFactor(), 0.0))
			continue;
		// Get the V size
		firstSize = vReactant.getSize();
		// Consider product with every HeV cluster.
		for (auto const& heVMapItem : getAll(ReactantType::PSIMixed)) {

			auto& heVReactant = *(heVMapItem.second);

			// Get its composition
			auto& comp = heVReactant.getComposition();
			// Create the composition of the potential product
			int newNumHe = comp[toCompIdx(Species::He)];
			int newNumD = comp[toCompIdx(Species::D)];
			int newNumT = comp[toCompIdx(Species::T)];
			int newNumV = comp[toCompIdx(Species::V)] + firstSize;

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::D)] = newNumD;
			newComp[toCompIdx(Species::T)] = newNumT;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::PSIMixed, newComp);

			// Check if the product can be a super cluster
			if (!product) {
				product = getSuperFromComp(newNumHe, newNumD, newNumT, newNumV);
			}
			// Check that the reaction can occur
			if (product
					&& (vReactant.getDiffusionFactor() > 0.0
							|| heVReactant.getDiffusionFactor() > 0.0)) {

				int a[4] = { newNumHe, newNumD, newNumT, newNumV };
				defineProductionReaction(vReactant, heVReactant, *product, a);
			}
		}

		// Consider product with each super cluster
		for (auto const& superMapItem : getAll(ReactantType::PSISuper)) {
			auto& superCluster =
					static_cast<PSISuperCluster&>(*(superMapItem.second));
			// Loop on the potential products
			for (auto const& superMapItemProd : getAll(ReactantType::PSISuper)) {
				auto& superProd =
						static_cast<PSISuperCluster&>(*(superMapItemProd.second));

				// Check if the reactions overlap
				auto& vCluster = static_cast<PSICluster&>(vReactant);
				bool overlap = checkOverlap(vCluster, superCluster, superProd);

				// Skip if they don't overlap
				if (!overlap)
					continue;

				// Check if the super clusters are full
				if (superCluster.isFull() && superProd.isFull()) {
					// This method will check if the reaction is possible and then add it to the list
					defineAnaProductionReactions(vReactant, superCluster,
							superProd);
				} else {
					std::vector<PendingProductionReactionInfo> prInfos;

					// Get the coordinates of the reactant and loop on them
					auto coords = superCluster.getCoordList();
					for (auto const& pair : coords) {
						// Assume the product can only be a super cluster here
						int newNumHe = std::get<0>(pair);
						int newNumD = std::get<1>(pair);
						int newNumT = std::get<2>(pair);
						int newNumV = std::get<3>(pair) + firstSize;
						if (superProd.isIn(newNumHe, newNumD, newNumT, newNumV)
								&& (vReactant.getDiffusionFactor() > 0.0
										|| superCluster.getDiffusionFactor()
												> 0.0)) {
							// Note that current reactant reacts with
							// current superCluster to produce product,
							// according to current parameters.
							int a[4] = { newNumHe, newNumD, newNumT, newNumV };
							int b[4] = { std::get<0>(pair), std::get<1>(pair),
									std::get<2>(pair), std::get<3>(pair) };
							prInfos.emplace_back(superProd, a, b);
						}
					}
					// Now that we know how current reactant reacts with
					// current superCluster, create the production
					// reaction(s) for them.
					if (prInfos.size() > 0)
						defineProductionReactions(vReactant, superCluster,
								prInfos);
				}
			}
		}
	}

	// Helium-Vacancy clustering
	// He_a + V_b --> (He_a)(V_b)
	// Consider each He cluster.
	for (auto const& heMapItem : getAll(ReactantType::He)) {

		auto& heReactant = *(heMapItem.second);

		// Get its size
		firstSize = heReactant.getSize();
		// Consider product with each V cluster.
		for (auto const& vMapItem : getAll(ReactantType::V)) {

			auto& vReactant = *(vMapItem.second);

			// Get its size
			secondSize = vReactant.getSize();
			// Create the composition of the potential product
			int newNumHe = firstSize;
			int newNumV = secondSize;

			// Get the product
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::PSIMixed, newComp);

			// Check if the product can be a super cluster
			if (!product) {
				product = getSuperFromComp(newNumHe, 0, 0, newNumV);
			}
			// Check that the reaction can occur
			if (product
					&& (heReactant.getDiffusionFactor() > 0.0
							|| vReactant.getDiffusionFactor() > 0.0)) {

				int a[4] = { newNumHe, 0, 0, newNumV };
				defineProductionReaction(heReactant, vReactant, *product, a);
			}
		}
	}

	// Vacancy reduction by Interstitial absorption in Mixed clusters
	// (He_a)()(V_b) + (I_c) --> (He_a)()[V_(b-c)]
	// Consider each I cluster
	for (auto const& iMapItem : getAll(ReactantType::I)) {

		auto& iReactant = *(iMapItem.second);
		// Skip if it can't diffuse
		if (xolotlCore::equal(iReactant.getDiffusionFactor(), 0.0))
			continue;

		// Get its size
		firstSize = iReactant.getSize();
		// Consider product with each HeV cluster.
		for (auto const& heVMapItem : getAll(ReactantType::PSIMixed)) {

			auto& heVReactant = *(heVMapItem.second);

			// Get its composition
			auto& comp = heVReactant.getComposition();
			// The product can be He or HeV
			IReactant * product = nullptr;
			if (comp[toCompIdx(Species::V)] == firstSize) {
				// The product can be He
				if (comp[toCompIdx(Species::D)] == 0
						&& comp[toCompIdx(Species::T)] == 0)
					product = get(Species::He, comp[toCompIdx(Species::He)]);
				// The product can be D
				if (comp[toCompIdx(Species::He)] == 0
						&& comp[toCompIdx(Species::T)] == 0)
					product = get(Species::D, comp[toCompIdx(Species::D)]);
				// The product can be T
				if (comp[toCompIdx(Species::D)] == 0
						&& comp[toCompIdx(Species::He)] == 0)
					product = get(Species::T, comp[toCompIdx(Species::T)]);
			} else {
				// The product is HeV
				// Create the composition of the potential product
				IReactant::Composition newComp;
				newComp[toCompIdx(Species::He)] = comp[toCompIdx(Species::He)];
				newComp[toCompIdx(Species::D)] = comp[toCompIdx(Species::D)];
				newComp[toCompIdx(Species::T)] = comp[toCompIdx(Species::T)];
				newComp[toCompIdx(Species::V)] = comp[toCompIdx(Species::V)]
						- firstSize;
				// Get the product
				product = get(ReactantType::PSIMixed, newComp);
			}
			// Check that the reaction can occur
			if (product
					&& (iReactant.getDiffusionFactor() > 0.0
							|| heVReactant.getDiffusionFactor() > 0.0)) {

				defineProductionReaction(iReactant, heVReactant, *product);
			}
		}

		// Consider product with all super clusters.
		for (auto const& superMapItem : getAll(ReactantType::PSISuper)) {
			auto& superCluster =
					static_cast<PSISuperCluster&>(*(superMapItem.second));
			std::vector<PendingProductionReactionInfo> prInfos;
			// Loop on the potential super products
			for (auto const& superMapItemProd : getAll(ReactantType::PSISuper)) {
				auto& superProd =
						static_cast<PSISuperCluster&>(*(superMapItemProd.second));

				// Check if the reactions overlap
				auto& iCluster = static_cast<PSICluster&>(iReactant);
				bool overlap = checkOverlap(iCluster, superCluster, superProd);

				// Skip if they don't overlap
				if (!overlap)
					continue;

				// Check if the super clusters are full
				if (superCluster.isFull() && superProd.isFull()) {
					// This method will check if the reaction is possible and then add it to the list
					defineAnaProductionReactions(iReactant, superCluster,
							superProd);
				} else {
					std::vector<PendingProductionReactionInfo> prInfos;

					// Get the coordinates of the reactant and loop on them
					auto coords = superCluster.getCoordList();
					for (auto const& pair : coords) {
						// Assume the product can only be a super cluster here
						int newNumHe = std::get<0>(pair);
						int newNumD = std::get<1>(pair);
						int newNumT = std::get<2>(pair);
						int newNumV = std::get<3>(pair) - firstSize;
						if (superProd.isIn(newNumHe, newNumD, newNumT, newNumV)
								&& (iReactant.getDiffusionFactor() > 0.0
										|| superCluster.getDiffusionFactor()
												> 0.0)) {
							// Note that current reactant reacts with
							// current superCluster to produce product,
							// according to current parameters.
							int a[4] = { newNumHe, newNumD, newNumT, newNumV };
							int b[4] = { std::get<0>(pair), std::get<1>(pair),
									std::get<2>(pair), std::get<3>(pair) };
							prInfos.emplace_back(superProd, a, b);
						}
					}
					// Now that we know how current reactant reacts with
					// current superCluster, create the production
					// reaction(s) for them.
					if (prInfos.size() > 0)
						defineProductionReactions(iReactant, superCluster,
								prInfos);
				}
			}

			// Get the coordinates of the reactant and loop on them
			auto coords = superCluster.getCoordList();
			for (auto const& pair : coords) {
				// The product might be mixed or He or D or T
				int newNumHe = std::get<0>(pair);
				int newNumD = std::get<1>(pair);
				int newNumT = std::get<2>(pair);
				int newNumV = std::get<3>(pair) - firstSize;

				// Get the product
				IReactant* product = nullptr;
				if (newNumV == 0) {
					// Check if it is a single product
					if ((newNumHe > 0) + (newNumD > 0) + (newNumT > 0) > 1) {
						// Nothing happens, no reaction
						continue;
					} else {
						if (newNumHe > 0)
							// The product is He
							product = get(Species::He, newNumHe);
						if (newNumD > 0)
							// The product is D
							product = get(Species::D, newNumD);
						if (newNumT > 0)
							// The product is T
							product = get(Species::T, newNumT);
					}
				} else {
					// Create the composition of the potential product
					IReactant::Composition newComp;
					newComp[toCompIdx(Species::He)] = newNumHe;
					newComp[toCompIdx(Species::D)] = newNumD;
					newComp[toCompIdx(Species::T)] = newNumT;
					newComp[toCompIdx(Species::V)] = newNumV;
					product = get(ReactantType::PSIMixed, newComp);
					// Don't test super product because it was already taken care of
				}
				// Check that the reaction can occur
				if (product
						&& (iReactant.getDiffusionFactor() > 0.0
								|| superCluster.getDiffusionFactor() > 0.0)) {
					int a[4] = { newNumHe, newNumD, newNumT, newNumV };
					int b[4] = { std::get<0>(pair), std::get<1>(pair), std::get<
							2>(pair), std::get<3>(pair) };
					prInfos.emplace_back(*product, a, b);
				}
			}

			// Now that we know how current reactant interacts with
			// current supercluster, define its production reactions
			// according to given parameters.
			if (prInfos.size() > 0)
				defineProductionReactions(iReactant, superCluster, prInfos);

		}
	}

	// Helium clustering leading to trap mutation
	// He_a + He_b --> [He_(a+b)](V_c) + I_c
	// Loop on the He clusters
	for (auto firstIt = getAll(ReactantType::He).begin();
			firstIt != getAll(ReactantType::He).end(); firstIt++) {
		// Get the reactant
		auto& firstReactant = static_cast<PSICluster&>(*(firstIt->second));
		// Get its size
		firstSize = firstReactant.getSize();
		// Loop on the second He cluster starting at the same pointer to avoid double counting
		for (auto secondIt = firstIt;
				secondIt != getAll(ReactantType::He).end(); secondIt++) {
			// Get the reactant
			auto& secondReactant = static_cast<PSICluster&>(*(secondIt->second));
			// Get its size
			secondSize = secondReactant.getSize();
			// Get the simple product
			productSize = firstSize + secondSize;
			auto product = get(toSpecies(ReactantType::He), productSize);
			// Doesn't do anything if the product exist
			if (product)
				continue;

			// Trap mutation is happening
			// Loop on the possible I starting by the smallest
			for (auto iSize = 1; iSize <= maxI; iSize++) {
				auto iReactant = get(toSpecies(ReactantType::I), iSize);
				// Create the composition of the potential product
				IReactant::Composition newComp;
				int newNumHe = firstSize + secondSize;
				int newNumV = iSize;
				newComp[toCompIdx(Species::He)] = newNumHe;
				newComp[toCompIdx(Species::V)] = newNumV;
				product = get(ReactantType::PSIMixed, newComp);

				// Check if the product can be a super cluster
				if (!product) {
					// Check if it is a super cluster from the map
					product = getSuperFromComp(newNumHe, 0, 0, newNumV);
				}

				// Check that the reaction can occur
				if (product
						&& (firstReactant.getDiffusionFactor() > 0.0
								|| secondReactant.getDiffusionFactor() > 0.0)) {

					int a[4] = { newNumHe, 0, 0, newNumV };
					defineProductionReaction(firstReactant, secondReactant,
							*iReactant, a);
					// This is a reaction with two products so we need to tell the other product
					// it is participating too
					int b[4] = { };
					defineProductionReaction(firstReactant, secondReactant,
							*product, a, b, true);

					// Stop the loop on I clusters here
					break;
				}
			}
		}
	}

	// Helium absorption by Mixed leading to trap mutation
	// (He_a)()(V_b) + He_c --> [He_(a+c)]()[V_(b+d)] + I_d
	// Loop on the He clusters
	for (auto const& heMapItem : getAll(ReactantType::He)) {
		auto& heReactant = *(heMapItem.second);

		// Skip if it can't diffuse
		if (xolotlCore::equal(heReactant.getDiffusionFactor(), 0.0))
			continue;
		// Get its size
		firstSize = heReactant.getSize();

		// Loop on the HeV clusters
		for (auto const& heVMapItem : getAll(ReactantType::PSIMixed)) {

			auto& heVReactant = *(heVMapItem.second);

			// Get its composition
			auto& comp = heVReactant.getComposition();
			// Create the composition of the potential product
			int newNumHe = comp[toCompIdx(Species::He)] + firstSize;
			int newNumD = comp[toCompIdx(Species::D)];
			int newNumT = comp[toCompIdx(Species::T)];
			int newNumV = comp[toCompIdx(Species::V)];

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::D)] = newNumD;
			newComp[toCompIdx(Species::T)] = newNumT;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::PSIMixed, newComp);
			// Doesn't do anything if the product exist
			if (product)
				continue;

			// Trap mutation is happening
			// Loop on the possible I starting by the smallest
			for (auto iSize = 1; iSize <= maxI; iSize++) {
				auto iReactant = get(toSpecies(ReactantType::I), iSize);
				// Create the composition of the potential product
				newComp[toCompIdx(Species::V)] = newNumV + iSize;
				product = get(ReactantType::PSIMixed, newComp);

				// Check if the product can be a super cluster
				if (!product) {
					// Check if it is a super cluster from the map
					product = getSuperFromComp(newNumHe, newNumD, newNumT,
							newNumV + iSize);
				}

				// Check that the reaction can occur
				if (product
						&& (heReactant.getDiffusionFactor() > 0.0
								|| heVReactant.getDiffusionFactor() > 0.0)) {

					int a[4] = { newNumHe, newNumD, newNumT, newNumV + iSize };
					defineProductionReaction(heReactant, heVReactant,
							*iReactant, a);
					// This is a reaction with two products so we need to tell the other product
					// it is participating too
					int b[4] = { };
					defineProductionReaction(heReactant, heVReactant, *product,
							a, b, true);

					// Stop the loop on I clusters here
					break;
				}
			}
		}

		// Consider product with all super clusters.
		for (auto const& superMapItem : getAll(ReactantType::PSISuper)) {
			auto& superCluster =
					static_cast<PSISuperCluster&>(*(superMapItem.second));
			std::vector<PendingProductionReactionInfo> prInfos1;
			std::vector<PendingProductionReactionInfo> prInfos2;

			// Get the coordinates of the reactant and loop on them
			auto coords = superCluster.getCoordList();
			for (auto const& pair : coords) {
				// The product might be mixed or He or D or T
				int newNumHe = std::get<0>(pair) + firstSize;
				int newNumD = std::get<1>(pair);
				int newNumT = std::get<2>(pair);
				int newNumV = std::get<3>(pair);

				// Get the product
				IReactant* product = getSuperFromComp(newNumHe, newNumD,
						newNumT, newNumV);
				// Skip if the product exists because we want trap mutation
				if (product)
					continue;

				// Trap mutation is happening
				// Loop on the possible I starting by the smallest
				for (auto iSize = 1; iSize <= maxI; iSize++) {
					auto iReactant = get(toSpecies(ReactantType::I), iSize);
					// Update the composition of the potential product
					product = getSuperFromComp(newNumHe, newNumD, newNumT,
							newNumV + iSize);

					// Check that the reaction can occur
					if (product && heReactant.getDiffusionFactor() > 0.0) {
						int a[4] =
								{ newNumHe, newNumD, newNumT, newNumV + iSize };
						int b[4] = { std::get<0>(pair), std::get<1>(pair),
								std::get<2>(pair), std::get<3>(pair) };
						prInfos1.emplace_back(*product, a, b);
						a[0] = 0, a[1] = 0, a[2] = 0, a[3] = 0;
						prInfos2.emplace_back(*iReactant, a, b);

						// Stop the loop on I clusters here
						break;
					}
				}
			}

			// Now that we know how current reactant interacts with
			// current supercluster, define its production reactions
			// according to given parameters.
			if (prInfos1.size() > 0 || prInfos2.size() > 0) {
				defineProductionReactions(heReactant, superCluster, prInfos2);
				defineProductionReactions(heReactant, superCluster, prInfos1,
						true);
			}
		}
	}

	// Vacancy-Interstitial annihilation
	// I_a + V_b
	//        --> I_(a-b), if a > b
	//        --> V_(b-a), if a < b
	//        --> 0, if a = b
	// Consider all I clusters.
	for (auto const& iMapItem : getAll(ReactantType::I)) {

		auto& iReactant = *(iMapItem.second);

		// Get its size
		firstSize = iReactant.getSize();

		// Consider product with each V cluster.
		for (auto const& vMapItem : getAll(ReactantType::V)) {

			auto& vReactant = *(vMapItem.second);

			// Get its size
			secondSize = vReactant.getSize();
			// Check the possibilities
			if (firstSize > secondSize) {
				// Get the product
				productSize = firstSize - secondSize;
				auto product = get(Species::I, productSize);
				// Check that the reaction can occur
				if (product
						&& (iReactant.getDiffusionFactor() > 0.0
								|| vReactant.getDiffusionFactor() > 0.0)) {

					defineAnnihilationReaction(iReactant, vReactant, *product);
				}
			} else if (firstSize < secondSize) {
				// Get the product
				productSize = secondSize - firstSize;
				auto product = get(Species::V, productSize);
				// Check that the reaction can occur
				if (product
						&& (iReactant.getDiffusionFactor() > 0.0
								|| vReactant.getDiffusionFactor() > 0.0)) {

					defineAnnihilationReaction(iReactant, vReactant, *product);
				}

			} else {
				// Annihilation
				// Check that the reaction can occur
				if ((iReactant.getDiffusionFactor() > 0.0
						|| vReactant.getDiffusionFactor() > 0.0)) {

					defineCompleteAnnihilationReaction(iReactant, vReactant);
				}
			}
		}
	}

	// Deuterium absorption by Mixed clusters
	// D_(a) + (D_b)()(V_c) --> [D_(a+b)]()(V_c)
	// Consider each D reactant.
	for (auto const& dMapItem : getAll(ReactantType::D)) {

		auto& dReactant = *(dMapItem.second);

		// Skip if it can't diffuse
		if (xolotlCore::equal(dReactant.getDiffusionFactor(), 0.0))
			continue;
		// Get its size
		firstSize = dReactant.getSize();

		// Consider product with each HeV cluster
		for (auto const& heVMapItem : getAll(ReactantType::PSIMixed)) {

			auto& heVReactant = *(heVMapItem.second);

			// Get its composition
			auto& comp = heVReactant.getComposition();
			// Create the composition of the potential product
			int newNumHe = comp[toCompIdx(Species::He)];
			int newNumD = comp[toCompIdx(Species::D)] + firstSize;
			int newNumT = comp[toCompIdx(Species::T)];
			int newNumV = comp[toCompIdx(Species::V)];

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::D)] = newNumD;
			newComp[toCompIdx(Species::T)] = newNumT;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::PSIMixed, newComp);

			// Check that the reaction can occur
			if (product
					&& (dReactant.getDiffusionFactor() > 0.0
							|| heVReactant.getDiffusionFactor() > 0.0)) {

				int a[4] = { newNumHe, newNumD, newNumT, newNumV };
				defineProductionReaction(dReactant, heVReactant, *product, a);
			}
		}

		// Consider product with each super cluster
		for (auto const& superMapItem : getAll(ReactantType::PSISuper)) {
			auto& superCluster =
					static_cast<PSISuperCluster&>(*(superMapItem.second));
			// Loop on the potential products
			for (auto const& superMapItemProd : getAll(ReactantType::PSISuper)) {
				auto& superProd =
						static_cast<PSISuperCluster&>(*(superMapItemProd.second));

				// Check if the reactions overlap
				auto& dCluster = static_cast<PSICluster&>(dReactant);
				bool overlap = checkOverlap(dCluster, superCluster, superProd);

				// Skip if they don't overlap
				if (!overlap)
					continue;

				// Check if the super clusters are full
				if (superCluster.isFull() && superProd.isFull()) {
					// This method will check if the reaction is possible and then add it to the list
					defineAnaProductionReactions(dReactant, superCluster,
							superProd);
				} else {
					std::vector<PendingProductionReactionInfo> prInfos;

					// Get the coordinates of the reactant and loop on them
					auto coords = superCluster.getCoordList();
					for (auto const& pair : coords) {
						// Assume the product can only be a super cluster here
						int newNumHe = std::get<0>(pair);
						int newNumD = std::get<1>(pair) + firstSize;
						int newNumT = std::get<2>(pair);
						int newNumV = std::get<3>(pair);
						if (superProd.isIn(newNumHe, newNumD, newNumT, newNumV)
								&& (dReactant.getDiffusionFactor() > 0.0
										|| superCluster.getDiffusionFactor()
												> 0.0)) {
							// Note that current reactant reacts with
							// current superCluster to produce product,
							// according to current parameters.
							int a[4] = { newNumHe, newNumD, newNumT, newNumV };
							int b[4] = { std::get<0>(pair), std::get<1>(pair),
									std::get<2>(pair), std::get<3>(pair) };
							prInfos.emplace_back(superProd, a, b);
						}
					}
					// Now that we know how current reactant reacts with
					// current superCluster, create the production
					// reaction(s) for them.
					if (prInfos.size() > 0)
						defineProductionReactions(dReactant, superCluster,
								prInfos);
				}
			}
		}
	}

	// Tritium absorption by Mixed clusters
	// T_(a) + (T_b)()(V_c) --> [T_(a+b)]()(V_c)
	// Consider each T reactant.
	for (auto const& tMapItem : getAll(ReactantType::T)) {

		auto& tReactant = *(tMapItem.second);

		// Skip if it can't diffuse
		if (xolotlCore::equal(tReactant.getDiffusionFactor(), 0.0))
			continue;
		// Get its size
		firstSize = tReactant.getSize();

		// Consider product with each Mixed cluster
		for (auto const& heVMapItem : getAll(ReactantType::PSIMixed)) {

			auto& heVReactant = *(heVMapItem.second);

			// Get its composition
			auto& comp = heVReactant.getComposition();
			// Create the composition of the potential product
			int newNumHe = comp[toCompIdx(Species::He)];
			int newNumD = comp[toCompIdx(Species::D)];
			int newNumT = comp[toCompIdx(Species::T)] + firstSize;
			int newNumV = comp[toCompIdx(Species::V)];

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::D)] = newNumD;
			newComp[toCompIdx(Species::T)] = newNumT;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::PSIMixed, newComp);

			// Check that the reaction can occur
			if (product
					&& (tReactant.getDiffusionFactor() > 0.0
							|| heVReactant.getDiffusionFactor() > 0.0)) {

				int a[4] = { newNumHe, newNumD, newNumT, newNumV };
				defineProductionReaction(tReactant, heVReactant, *product, a);
			}
		}

		// Consider product with each super cluster
		for (auto const& superMapItem : getAll(ReactantType::PSISuper)) {
			auto& superCluster =
					static_cast<PSISuperCluster&>(*(superMapItem.second));
			// Loop on the potential products
			for (auto const& superMapItemProd : getAll(ReactantType::PSISuper)) {
				auto& superProd =
						static_cast<PSISuperCluster&>(*(superMapItemProd.second));

				// Check if the reactions overlap
				auto& tCluster = static_cast<PSICluster&>(tReactant);
				bool overlap = checkOverlap(tCluster, superCluster, superProd);

				// Skip if they don't overlap
				if (!overlap)
					continue;

				// Check if the super clusters are full
				if (superCluster.isFull() && superProd.isFull()) {
					// This method will check if the reaction is possible and then add it to the list
					defineAnaProductionReactions(tReactant, superCluster,
							superProd);
				} else {
					std::vector<PendingProductionReactionInfo> prInfos;

					// Get the coordinates of the reactant and loop on them
					auto coords = superCluster.getCoordList();
					for (auto const& pair : coords) {
						// Assume the product can only be a super cluster here
						int newNumHe = std::get<0>(pair);
						int newNumD = std::get<1>(pair);
						int newNumT = std::get<2>(pair) + firstSize;
						int newNumV = std::get<3>(pair);
						if (superProd.isIn(newNumHe, newNumD, newNumT, newNumV)
								&& (tReactant.getDiffusionFactor() > 0.0
										|| superCluster.getDiffusionFactor()
												> 0.0)) {
							// Note that current reactant reacts with
							// current superCluster to produce product,
							// according to current parameters.
							int a[4] = { newNumHe, newNumD, newNumT, newNumV };
							int b[4] = { std::get<0>(pair), std::get<1>(pair),
									std::get<2>(pair), std::get<3>(pair) };
							prInfos.emplace_back(superProd, a, b);
						}
					}
					// Now that we know how current reactant reacts with
					// current superCluster, create the production
					// reaction(s) for them.
					if (prInfos.size() > 0)
						defineProductionReactions(tReactant, superCluster,
								prInfos);
				}
			}
		}
	}

	// Deuterium-Vacancy clustering
	// D_a + V_b --> (D_a)(V_b)
	// Consider each D cluster.
	for (auto const& dMapItem : getAll(ReactantType::D)) {

		auto& dReactant = *(dMapItem.second);

		// Get its size
		firstSize = dReactant.getSize();
		// Consider product with each V cluster.
		for (auto const& vMapItem : getAll(ReactantType::V)) {

			auto& vReactant = *(vMapItem.second);

			// Get its size
			secondSize = vReactant.getSize();
			// Create the composition of the potential product
			int newNumD = firstSize;
			int newNumV = secondSize;

			// Get the product
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::D)] = newNumD;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::PSIMixed, newComp);

			// Check if the product can be a super cluster
			if (!product) {
				product = getSuperFromComp(0, newNumD, 0, newNumV);
			}

			// Check that the reaction can occur
			if (product
					&& (dReactant.getDiffusionFactor() > 0.0
							|| vReactant.getDiffusionFactor() > 0.0)) {

				int a[4] = { 0, newNumD, 0, newNumV };
				defineProductionReaction(dReactant, vReactant, *product, a);
			}
		}
	}

	// Tritium-Vacancy clustering
	// T_a + V_b --> (T_a)(V_b)
	// Consider each T cluster.
	for (auto const& tMapItem : getAll(ReactantType::T)) {

		auto& tReactant = *(tMapItem.second);

		// Get its size
		firstSize = tReactant.getSize();
		// Consider product with each V cluster.
		for (auto const& vMapItem : getAll(ReactantType::V)) {

			auto& vReactant = *(vMapItem.second);

			// Get its size
			secondSize = vReactant.getSize();
			// Create the composition of the potential product
			int newNumT = firstSize;
			int newNumV = secondSize;

			// Get the product
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::T)] = newNumT;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::PSIMixed, newComp);

			// Check if the product can be a super cluster
			if (!product) {
				product = getSuperFromComp(0, 0, newNumT, newNumV);
			}

			// Check that the reaction can occur
			if (product
					&& (tReactant.getDiffusionFactor() > 0.0
							|| vReactant.getDiffusionFactor() > 0.0)) {

				int a[4] = { 0, 0, newNumT, newNumV };
				defineProductionReaction(tReactant, vReactant, *product, a);
			}
		}
	}

	// Helium absorption by HeI clusters
	// He_(a) + (He_b)(I_c) --> [He_(a+b)](I_c)
	// Consider each He reactant.
	for (auto const& heMapItem : getAll(ReactantType::He)) {

		auto& heReactant = *(heMapItem.second);

		// Skip if it can't diffuse
		if (xolotlCore::equal(heReactant.getDiffusionFactor(), 0.0))
			continue;
		// Get its size
		firstSize = heReactant.getSize();

		// Consider product with each HeI cluster
		for (auto const& heIMapItem : getAll(ReactantType::HeI)) {

			auto& heIReactant = *(heIMapItem.second);

			// Get its composition
			auto& comp = heIReactant.getComposition();
			// Create the composition of the potential product
			int newNumHe = comp[toCompIdx(Species::He)] + firstSize;
			int newNumI = comp[toCompIdx(Species::I)];

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::I)] = newNumI;
			auto product = get(ReactantType::HeI, newComp);

			// Check that the reaction can occur
			if (product
					&& (heReactant.getDiffusionFactor() > 0.0
							|| heIReactant.getDiffusionFactor() > 0.0)) {

				int a[4] = { newNumHe, 0, 0, newNumI };
				defineProductionReaction(heReactant, heIReactant, *product, a);
			}
		}
	}

	// Interstitial absorption by HeI clusters
	// (He_a)(I_b) + I_c --> (He_a)[I_(b+c)]
	// Consider each I cluster.
	for (auto const& iMapItem : getAll(ReactantType::I)) {

		auto& iReactant = *(iMapItem.second);

		// Skip if it can't diffuse
		if (xolotlCore::equal(iReactant.getDiffusionFactor(), 0.0))
			continue;
		// Get the I size
		firstSize = iReactant.getSize();
		// Consider product with every HeI cluster.
		for (auto const& heIMapItem : getAll(ReactantType::HeI)) {

			auto& heIReactant = *(heIMapItem.second);

			// Get its composition
			auto& comp = heIReactant.getComposition();
			// Create the composition of the potential product
			int newNumHe = comp[toCompIdx(Species::He)];
			int newNumI = comp[toCompIdx(Species::I)] + firstSize;

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::I)] = newNumI;
			auto product = get(ReactantType::HeI, newComp);

			// Check that the reaction can occur
			if (product
					&& (iReactant.getDiffusionFactor() > 0.0
							|| heIReactant.getDiffusionFactor() > 0.0)) {

				int a[4] = { newNumHe, 0, 0, newNumI };
				defineProductionReaction(iReactant, heIReactant, *product, a);
			}
		}
	}

	// Helium-Interstitial clustering
	// He_a + I_b --> (He_a)(I_b)
	// Consider each He cluster.
	for (auto const& heMapItem : getAll(ReactantType::He)) {

		auto& heReactant = *(heMapItem.second);

		// Get its size
		firstSize = heReactant.getSize();
		// Consider product with each I cluster.
		for (auto const& iMapItem : getAll(ReactantType::I)) {

			auto& iReactant = *(iMapItem.second);

			// Get its size
			secondSize = iReactant.getSize();
			// Create the composition of the potential product
			int newNumHe = firstSize;
			int newNumI = secondSize;

			// Get the product
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::I)] = newNumI;
			auto product = get(ReactantType::HeI, newComp);

			// Check that the reaction can occur
			if (product
					&& (heReactant.getDiffusionFactor() > 0.0
							|| iReactant.getDiffusionFactor() > 0.0)) {

				int a[4] = { newNumHe, 0, 0, newNumI };
				defineProductionReaction(heReactant, iReactant, *product, a);
			}
		}
	}

	// Interstitial reduction by vacancy absorption in HeI clusters
	// (He_a)(I_b) + (V_c) --> (He_a)[I_(b-c)]
	// Consider each V cluster
	for (auto const& vMapItem : getAll(ReactantType::V)) {

		auto& vReactant = *(vMapItem.second);

		// Get its size
		firstSize = vReactant.getSize();
		// Consider product with each HeI cluster.
		for (auto const& heIMapItem : getAll(ReactantType::HeI)) {

			auto& heIReactant = *(heIMapItem.second);

			// Get its composition
			auto& comp = heIReactant.getComposition();
			// The product can be He or HeV
			IReactant * product = nullptr;
			if (comp[toCompIdx(Species::I)] == firstSize) {
				// The product is He
				product = get(Species::He, comp[toCompIdx(Species::He)]);
			} else {
				// The product is HeI
				// Create the composition of the potential product
				IReactant::Composition newComp;
				newComp[toCompIdx(Species::He)] = comp[toCompIdx(Species::He)];
				newComp[toCompIdx(Species::I)] = comp[toCompIdx(Species::I)]
						- firstSize;
				// Get the product
				product = get(ReactantType::HeI, newComp);
			}
			// Check that the reaction can occur
			if (product
					&& (vReactant.getDiffusionFactor() > 0.0
							|| heIReactant.getDiffusionFactor() > 0.0)) {

				defineProductionReaction(vReactant, heIReactant, *product);
			}
		}
	}

	return;
}

bool PSIClusterReactionNetwork::canDissociate(IReactant& emittingReactant,
		ProductionReaction& reaction) const {
	// Assume reaction can dissociate by default.
	bool ret = true;

	// Check if at least one of the potentially emitted cluster is size one
	if (reaction.first.getSize() != 1 && reaction.second.getSize() != 1) {
		// Don't add the reverse reaction
		ret = false;
	}

	// Check for trap mutations
	if (emittingReactant.getType() == ReactantType::I
			&& reaction.first.getType() != ReactantType::I) {
		// Don't add the reverse reaction
		ret = false;
	}
	// (with XOR)
	if ((reaction.first.getType() == ReactantType::I)
			== !(reaction.second.getType() == ReactantType::I)) {
		// Don't add the reverse reaction
		ret = false;
	}

	return ret;
}

void PSIClusterReactionNetwork::checkForDissociation(
		IReactant& emittingReactant, ProductionReaction& reaction, int a[4],
		int b[4]) {

	// Check if reaction can dissociate.
	if (canDissociate(emittingReactant, reaction)) {
		// The dissociation can occur, so create a reaction for it.
		defineDissociationReaction(reaction, emittingReactant, a, b);
	}

	return;
}

void PSIClusterReactionNetwork::setTemperature(double temp) {
	ReactionNetwork::setTemperature(temp);

	computeRateConstants();

	return;
}

void PSIClusterReactionNetwork::reinitializeNetwork() {

	// Reset the Ids
	// std::for_each is guaranteed to visit reactants in order for C++11.
	int id = 0;
	std::for_each(allReactants.begin(), allReactants.end(),
			[&id](IReactant& currReactant) {
				id++;
				currReactant.setId(id);
				currReactant.setMomentId(id, 0);
				currReactant.setMomentId(id, 1);
				currReactant.setMomentId(id, 2);
				currReactant.setMomentId(id, 3);
			});

	// Get all the super clusters and loop on them
	// Have to use allReactants again to be sure the ordering is the same across plateforms
	std::for_each(allReactants.begin(), allReactants.end(),
			[&id, this](IReactant& currReactant) {

				if (currReactant.getType() == ReactantType::PSISuper) {
					auto& currCluster = static_cast<PSISuperCluster&>(currReactant);
					// Loop on the axis
					for (int i = 1; i < psDim; i++) {
						id++;
						currCluster.setMomentId(id, indexList[i] - 1);
					}

					// Update the PSIMixed size
					IReactant::SizeType clusterSize = currCluster.getBounds(0).second
					+ currCluster.getBounds(3).second;
					if (clusterSize > maxClusterSizeMap[ReactantType::PSIMixed]) {
						maxClusterSizeMap[ReactantType::PSIMixed] = clusterSize;
					}
				}
			});

	return;
}

void PSIClusterReactionNetwork::reinitializeConnectivities() {

	// Reset connectivities of each reactant.
	std::for_each(allReactants.begin(), allReactants.end(),
			[](IReactant& currReactant) {
				currReactant.resetConnectivities();
			});

	return;
}

void PSIClusterReactionNetwork::updateConcentrationsFromArray(
		double * concentrations) {

	// Set the concentration on each reactant.
	std::for_each(allReactants.begin(), allReactants.end(),
			[&concentrations](IReactant& currReactant) {
				auto id = currReactant.getId() - 1;
				currReactant.setConcentration(concentrations[id]);
			});

	// Set the moments
	auto const& superTypeMap = getAll(ReactantType::PSISuper);
	std::for_each(superTypeMap.begin(), superTypeMap.end(),
			[&concentrations,this](const ReactantMap::value_type& currMapItem) {

				auto& cluster = static_cast<PSISuperCluster&>(*(currMapItem.second));

				cluster.setZerothMoment(concentrations[cluster.getId() - 1]);
				// Loop on the used moments
				for (int i = 1; i < psDim; i++) {
					cluster.setMoment(concentrations[cluster.getMomentId(indexList[i] - 1) - 1], indexList[i] - 1);
				}
			});

	return;
}

std::vector<std::vector<int> > PSIClusterReactionNetwork::getCompositionList() const {
	// Create the list that will be returned
	std::vector<std::vector<int> > compList;

	// Loop on all the reactants
	std::for_each(allReactants.begin(), allReactants.end(),
			[&compList](IReactant& currReactant) {
				// Get the composition
				auto comp = currReactant.getComposition();
				std::vector <int> compVec;
				compVec.push_back(comp[toCompIdx(Species::He)]);
				compVec.push_back(comp[toCompIdx(Species::D)]);
				compVec.push_back(comp[toCompIdx(Species::T)]);
				compVec.push_back(comp[toCompIdx(Species::V)]);
				compVec.push_back(comp[toCompIdx(Species::I)]);

				// Save the composition in the list
				compList.push_back(compVec);
			});

	return compList;
}

void PSIClusterReactionNetwork::getDiagonalFill(SparseFillMap& fillMap) {

	// Degrees of freedom is the total number of clusters in the network
	const int dof = getDOF();

	// Get the connectivity for each reactant
	std::for_each(allReactants.begin(), allReactants.end(),
			[&fillMap,&dof,this](const IReactant& reactant) {

				// Get the reactant's connectivity
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
					if(connectivity[j] == 1) {
						// TODO are fillMap and dFillmap the same?
						fillMap[id].emplace_back(j);
						columnIds.emplace_back(j);
					}
				}
				// Update the map
				dFillMap[id] = columnIds;
			});

	// Get the connectivity for each moment
	for (auto const& superMapItem : getAll(ReactantType::PSISuper)) {

		// Get the reactant and its connectivity
		auto const& reactant =
				static_cast<PSISuperCluster&>(*(superMapItem.second));
		auto const& connectivity = reactant.getConnectivity();
		auto connectivityLength = connectivity.size();

		// Loop on the axis
		for (int i = 1; i < psDim; i++) {
			// Get the helium moment id so that the connectivity can be lined up in
			// the proper column
			auto id = reactant.getMomentId(indexList[i] - 1) - 1;

			// Create the vector that will be inserted into the dFill map
			std::vector<int> columnIds;
			// Add it to the diagonal fill block
			for (int j = 0; j < connectivityLength; j++) {
				// Add a column id if the connectivity is equal to 1.
				if (connectivity[j] == 1) {
					fillMap[id].emplace_back(j);
					columnIds.emplace_back(j);
				}
				// Update the map
				dFillMap[id] = columnIds;
			}
		}
	}

	// Now that the dFillMap has been built, build inverse maps for each item.
	for (const auto& dFillMapItem : dFillMap) {
		auto rid = dFillMapItem.first;
		dFillInvMap[rid] = PartialsIdxMap();

		auto const& colIds = dFillMapItem.second;
		for (auto j = 0; j < colIds.size(); ++j) {
			dFillInvMap[rid][colIds[j]] = j;
		}
	}

	return;
}

double PSIClusterReactionNetwork::getTotalAtomConcentration(int i) {
	// Initial declarations
	double atomConc = 0.0;
	ReactantType type;

	// Switch on the index
	switch (i) {
	case 0:
		type = ReactantType::He;
		break;
	case 1:
		type = ReactantType::D;
		break;
	case 2:
		type = ReactantType::T;
		break;
	default:
		throw std::string("\nType not defined for getTotalAtomConcentration()");
		break;
	}

	// Sum over all He clusters.
	for (auto const& currMapItem : getAll(type)) {

		// Get the cluster and its composition
		auto const& cluster = *(currMapItem.second);
		double size = cluster.getSize();

		// Add the concentration times the He content to the total helium concentration
		atomConc += cluster.getConcentration() * size;
	}

	// Sum over all Mixed clusters.
	for (auto const& currMapItem : getAll(ReactantType::PSIMixed)) {

		// Get the cluster and its composition
		auto const& cluster = *(currMapItem.second);
		auto& comp = cluster.getComposition();

		// Add the concentration times the He content to the total helium concentration
		atomConc += cluster.getConcentration()
				* comp[toCompIdx(toSpecies(type))];
	}

	// Sum over all super clusters.
	for (auto const& currMapItem : getAll(ReactantType::PSISuper)) {

		// Get the cluster
		auto const& cluster =
				static_cast<PSISuperCluster&>(*(currMapItem.second));

		// Add its total atom concentration
		atomConc += cluster.getTotalAtomConcentration(i);
	}

	return atomConc;
}

double PSIClusterReactionNetwork::getTotalTrappedAtomConcentration(int i) {
	// Initial declarations
	double atomConc = 0.0;
	ReactantType type;

	// Switch on the index
	switch (i) {
	case 0:
		type = ReactantType::He;
		break;
	case 1:
		type = ReactantType::D;
		break;
	case 2:
		type = ReactantType::T;
		break;
	default:
		throw std::string("\nType not defined for getTotalAtomConcentration()");
		break;
	}

	// Sum over all Mixed clusters.
	for (auto const& currMapItem : getAll(ReactantType::PSIMixed)) {

		// Get the cluster and its composition
		auto const& cluster = *(currMapItem.second);
		auto& comp = cluster.getComposition();

		// Add the concentration times the He content to the total helium concentration
		atomConc += cluster.getConcentration()
				* comp[toCompIdx(toSpecies(type))];
	}

	// Sum over all super clusters.
	for (auto const& currMapItem : getAll(ReactantType::PSISuper)) {

		// Get the cluster
		auto const& cluster =
				static_cast<PSISuperCluster&>(*(currMapItem.second));

		// Add its total helium concentration helium concentration
		atomConc += cluster.getTotalAtomConcentration(i);
	}

	return atomConc;
}

double PSIClusterReactionNetwork::getTotalVConcentration() {
	// Initial declarations
	double vConc = 0.0;

	// Sum over all V clusters.
	for (auto const& currMapItem : getAll(ReactantType::V)) {
		// Get the cluster and its composition
		auto const& cluster = *(currMapItem.second);
		double size = cluster.getSize();

		// Add the concentration times the V content to the total vacancy concentration
		vConc += cluster.getConcentration() * size;
	}

	// Sum over all HeV clusters
	for (auto const& currMapItem : getAll(ReactantType::PSIMixed)) {
		// Get the cluster and its composition
		auto const& cluster = *(currMapItem.second);
		auto& comp = cluster.getComposition();

		// Add the concentration times the V content to the total vacancy concentration
		vConc += cluster.getConcentration() * comp[toCompIdx(Species::V)];
	}

	// Sum over all super clusters
	for (auto const& currMapItem : getAll(ReactantType::PSISuper)) {
		// Get the cluster
		auto const& cluster =
				static_cast<PSISuperCluster&>(*(currMapItem.second));

		// Add its total vacancy concentration
		vConc += cluster.getTotalVacancyConcentration();
	}

	return vConc;
}

double PSIClusterReactionNetwork::getTotalIConcentration() {
	// Initial declarations
	double iConc = 0.0;

	// Sum over all I clusters
	for (auto const& currMapItem : getAll(ReactantType::I)) {
		// Get the cluster and its composition
		auto const& cluster = *(currMapItem.second);
		double size = cluster.getSize();

		// Add the concentration times the I content to the total interstitial concentration
		iConc += cluster.getConcentration() * size;
	}

	return iConc;
}

void PSIClusterReactionNetwork::computeRateConstants() {
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

void PSIClusterReactionNetwork::computeAllFluxes(double *updatedConcOffset) {

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
	for (auto const& currMapItem : getAll(ReactantType::PSISuper)) {

		auto const& superCluster =
				static_cast<PSISuperCluster&>(*(currMapItem.second));

		// Loop on the axis
		for (int i = 1; i < psDim; i++) {

			// Compute the moment flux
			auto flux = superCluster.getMomentFlux(indexList[i] - 1);
			// Update the concentration of the cluster
			auto reactantIndex = superCluster.getMomentId(indexList[i] - 1) - 1;
			updatedConcOffset[reactantIndex] += flux;
		}
	}

	return;
}

void PSIClusterReactionNetwork::computeAllPartials(
		const std::vector<size_t>& startingIdx, const std::vector<int>& indices,
		std::vector<double>& vals) const {

	// Because we accumulate partials and we don't know which
	// of our reactants will be first to assign a value, we must start with
	// all partials values at zero.
	std::fill(vals.begin(), vals.end(), 0.0);

	// Initial declarations
	std::vector<double> clusterPartials(getDOF(), 0.0);

	// Make a vector of types for the non super clusters
	std::vector<ReactantType> typeVec { ReactantType::He, ReactantType::D,
			ReactantType::T, ReactantType::V, ReactantType::I,
			ReactantType::PSIMixed };
	// Loop on it
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {

		auto currType = *tvIter;

		// Consider all reactants of the current type.
		auto const& currTypeReactantMap = getAll(currType);

		// Update the column in the Jacobian that represents each normal reactant
		for (auto const& currMapItem : currTypeReactantMap) {

			auto const& reactant =
					static_cast<PSICluster&>(*(currMapItem.second));

			// Get the reactant index
			auto reactantIndex = reactant.getId() - 1;

			// Get the partial derivatives
			reactant.getPartialDerivatives(clusterPartials);
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

	// Update the column in the Jacobian that represents the moment for the super clusters
	// Create the partials container that is going to be used
	double** partials = new double*[5];
	auto const& superClusters = getAll(ReactantType::PSISuper);
	for (auto const& currMapItem : superClusters) {

		auto const& reactant =
				static_cast<PSISuperCluster&>(*(currMapItem.second));

		// Determine cluster's index into the size/indices/vals arrays.
		int reactantIndices[5] = { };
		reactantIndices[0] = reactant.getId() - 1;
		// Loop on the axis for the moments
		for (int i = 1; i < psDim; i++) {
			// Get the moment index
			reactantIndices[i] = reactant.getMomentId(indexList[i] - 1) - 1;
		}

		// Get the inverse mappings from dense DOF space to
		// the indices/vals arrays.
		PartialsIdxMap partialsIdxMap[5];
		for (int i = 0; i < psDim; i++) {
			partialsIdxMap[i] = dFillInvMap.at(reactantIndices[i]);
			partials[i] = &(vals[startingIdx[reactantIndices[i]]]);
		}

		// Have reactant compute its partial derivatives
		// to its correct locations within the vals array.
		reactant.computePartialDerivatives(partials, partialsIdxMap);
	}

	// Clear memory
	delete[] partials;

	return;
}

double PSIClusterReactionNetwork::computeBindingEnergy(
		const DissociationReaction& reaction) const {
// for the dissociation A --> B + C we need A binding energy
// E_b(A) = E_f(B) + E_f(C) - E_f(A) where E_f is the formation energy
	double bindingEnergy = reaction.first.getFormationEnergy()
			+ reaction.second.getFormationEnergy()
			- reaction.dissociating.getFormationEnergy();

// hydrogen cases
	if (reaction.dissociating.getType() == ReactantType::PSIMixed
			&& (reaction.first.getType() == ReactantType::D
					|| reaction.first.getType() == ReactantType::T
					|| reaction.second.getType() == ReactantType::D
					|| reaction.second.getType() == ReactantType::T)) {
		auto comp = reaction.dissociating.getComposition();
		int heSize = comp[toCompIdx(Species::He)];
		int vSize = comp[toCompIdx(Species::V)];
		int hSize = comp[toCompIdx(Species::D)] + comp[toCompIdx(Species::T)];

		if (vSize == 1) {
			switch (heSize) {
			case 0:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.21;
					break;
				case 2:
					bindingEnergy = 1.17;
					break;
				case 3:
					bindingEnergy = 1.05;
					break;
				case 4:
					bindingEnergy = 0.93;
					break;
				case 5:
					bindingEnergy = 0.85;
					break;
				case 6:
					bindingEnergy = 0.60;
					break;
				default:
					break;
				}
				break;
			case 1:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.00;
					break;
				case 2:
					bindingEnergy = 0.95;
					break;
				case 3:
					bindingEnergy = 0.90;
					break;
				case 4:
					bindingEnergy = 0.88;
					break;
				case 5:
					bindingEnergy = 0.80;
					break;
				case 6:
					bindingEnergy = 0.60;
					break;
				default:
					break;
				}
				break;
			case 2:
				switch (hSize) {
				case 1:
					bindingEnergy = 0.96;
					break;
				case 2:
					bindingEnergy = 0.92;
					break;
				case 3:
					bindingEnergy = 0.85;
					break;
				case 4:
					bindingEnergy = 0.84;
					break;
				case 5:
					bindingEnergy = 0.83;
					break;
				case 6:
					bindingEnergy = 0.50;
					break;
				default:
					break;
				}
				break;
			case 3:
				switch (hSize) {
				case 1:
					bindingEnergy = 0.86;
					break;
				case 2:
					bindingEnergy = 0.81;
					break;
				case 3:
					bindingEnergy = 0.69;
					break;
				case 4:
					bindingEnergy = 0.64;
					break;
				case 5:
					bindingEnergy = 0.65;
					break;
				case 6:
					bindingEnergy = 0.50;
					break;
				default:
					break;
				}
				break;
			case 4:
				switch (hSize) {
				case 1:
					bindingEnergy = 0.83;
					break;
				case 2:
					bindingEnergy = 0.80;
					break;
				case 3:
					bindingEnergy = 0.65;
					break;
				case 4:
					bindingEnergy = 0.60;
					break;
				case 5:
					bindingEnergy = 0.60;
					break;
				case 6:
					bindingEnergy = 0.55;
					break;
				default:
					break;
				}
				break;
			case 5:
				switch (hSize) {
				case 1:
					bindingEnergy = 0.83;
					break;
				case 2:
					bindingEnergy = 0.80;
					break;
				case 3:
					bindingEnergy = 0.60;
					break;
				case 4:
					bindingEnergy = 0.50;
					break;
				case 5:
					bindingEnergy = 0.50;
					break;
				case 6:
					bindingEnergy = 0.50;
					break;
				default:
					break;
				}
				break;
			case 6:
				switch (hSize) {
				case 1:
					bindingEnergy = 0.80;
					break;
				case 2:
					bindingEnergy = 0.70;
					break;
				case 3:
					bindingEnergy = 0.60;
					break;
				case 4:
					bindingEnergy = 0.50;
					break;
				case 5:
					bindingEnergy = 0.50;
					break;
				case 6:
					bindingEnergy = 0.50;
					break;
				default:
					break;
				}
				break;
			case 7:
				switch (hSize) {
				case 1:
					bindingEnergy = 0.80;
					break;
				case 2:
					bindingEnergy = 0.75;
					break;
				case 3:
					bindingEnergy = 0.65;
					break;
				case 4:
					bindingEnergy = 0.55;
					break;
				case 5:
					bindingEnergy = 0.55;
					break;
				case 6:
					bindingEnergy = 0.45;
					break;
				default:
					break;
				}
				break;
			case 8:
				switch (hSize) {
				case 1:
					bindingEnergy = 0.80;
					break;
				case 2:
					bindingEnergy = 0.80;
					break;
				case 3:
					bindingEnergy = 0.70;
					break;
				case 4:
					bindingEnergy = 0.65;
					break;
				case 5:
					bindingEnergy = 0.60;
					break;
				case 6:
					bindingEnergy = 0.55;
					break;
				default:
					break;
				}
				break;
			case 9:
				switch (hSize) {
				case 1:
					bindingEnergy = 0.80;
					break;
				case 2:
					bindingEnergy = 0.80;
					break;
				case 3:
					bindingEnergy = 0.75;
					break;
				case 4:
					bindingEnergy = 0.70;
					break;
				case 5:
					bindingEnergy = 0.65;
					break;
				case 6:
					bindingEnergy = 0.60;
					break;
				default:
					break;
				}
				break;
			default:
				break;
			}
		} else if (vSize == 2) {
			switch (heSize) {
			case 0:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.63;
					break;
				case 2:
					bindingEnergy = 1.31;
					break;
				case 3:
					bindingEnergy = 1.25;
					break;
				case 4:
					bindingEnergy = 1.16;
					break;
				case 5:
					bindingEnergy = 1.00;
					break;
				case 6:
					bindingEnergy = 1.00;
					break;
				case 7:
					bindingEnergy = 0.95;
					break;
				case 8:
					bindingEnergy = 0.95;
					break;
				case 9:
					bindingEnergy = 0.75;
					break;
				case 10:
					bindingEnergy = 0.70;
					break;
				case 11:
					bindingEnergy = 0.65;
					break;
				default:
					break;
				}
				break;
			case 1:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.30;
					break;
				case 2:
					bindingEnergy = 1.30;
					break;
				case 3:
					bindingEnergy = 1.24;
					break;
				case 4:
					bindingEnergy = 1.08;
					break;
				case 5:
					bindingEnergy = 0.95;
					break;
				case 6:
					bindingEnergy = 0.95;
					break;
				case 7:
					bindingEnergy = 0.95;
					break;
				case 8:
					bindingEnergy = 0.95;
					break;
				case 9:
					bindingEnergy = 0.75;
					break;
				case 10:
					bindingEnergy = 0.70;
					break;
				case 11:
					bindingEnergy = 0.65;
					break;
				default:
					break;
				}
				break;
			case 2:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.15;
					break;
				case 2:
					bindingEnergy = 1.14;
					break;
				case 3:
					bindingEnergy = 1.11;
					break;
				case 4:
					bindingEnergy = 1.14;
					break;
				case 5:
					bindingEnergy = 0.95;
					break;
				case 6:
					bindingEnergy = 0.95;
					break;
				case 7:
					bindingEnergy = 0.95;
					break;
				case 8:
					bindingEnergy = 0.90;
					break;
				case 9:
					bindingEnergy = 0.75;
					break;
				case 10:
					bindingEnergy = 0.70;
					break;
				case 11:
					bindingEnergy = 0.65;
					break;
				default:
					break;
				}
				break;
			case 3:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.12;
					break;
				case 2:
					bindingEnergy = 1.06;
					break;
				case 3:
					bindingEnergy = 0.99;
					break;
				case 4:
					bindingEnergy = 0.99;
					break;
				case 5:
					bindingEnergy = 0.90;
					break;
				case 6:
					bindingEnergy = 0.95;
					break;
				case 7:
					bindingEnergy = 0.90;
					break;
				case 8:
					bindingEnergy = 0.90;
					break;
				case 9:
					bindingEnergy = 0.70;
					break;
				case 10:
					bindingEnergy = 0.70;
					break;
				case 11:
					bindingEnergy = 0.65;
					break;
				default:
					break;
				}
				break;
			case 4:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.10;
					break;
				case 2:
					bindingEnergy = 1.06;
					break;
				case 3:
					bindingEnergy = 0.99;
					break;
				case 4:
					bindingEnergy = 0.99;
					break;
				case 5:
					bindingEnergy = 0.90;
					break;
				case 6:
					bindingEnergy = 0.95;
					break;
				case 7:
					bindingEnergy = 0.90;
					break;
				case 8:
					bindingEnergy = 0.90;
					break;
				case 9:
					bindingEnergy = 0.70;
					break;
				case 10:
					bindingEnergy = 0.65;
					break;
				case 11:
					bindingEnergy = 0.65;
					break;
				default:
					break;
				}
				break;
			case 5:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.10;
					break;
				case 2:
					bindingEnergy = 1.05;
					break;
				case 3:
					bindingEnergy = 0.99;
					break;
				case 4:
					bindingEnergy = 0.99;
					break;
				case 5:
					bindingEnergy = 0.90;
					break;
				case 6:
					bindingEnergy = 0.90;
					break;
				case 7:
					bindingEnergy = 0.90;
					break;
				case 8:
					bindingEnergy = 0.90;
					break;
				case 9:
					bindingEnergy = 0.70;
					break;
				case 10:
					bindingEnergy = 0.65;
					break;
				case 11:
					bindingEnergy = 0.65;
					break;
				default:
					break;
				}
				break;
			case 6:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.10;
					break;
				case 2:
					bindingEnergy = 1.05;
					break;
				case 3:
					bindingEnergy = 0.99;
					break;
				case 4:
					bindingEnergy = 0.99;
					break;
				case 5:
					bindingEnergy = 0.90;
					break;
				case 6:
					bindingEnergy = 0.90;
					break;
				case 7:
					bindingEnergy = 0.90;
					break;
				case 8:
					bindingEnergy = 0.85;
					break;
				case 9:
					bindingEnergy = 0.70;
					break;
				case 10:
					bindingEnergy = 0.65;
					break;
				case 11:
					bindingEnergy = 0.60;
					break;
				default:
					break;
				}
				break;
			case 7:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.05;
					break;
				case 2:
					bindingEnergy = 1.00;
					break;
				case 3:
					bindingEnergy = 0.95;
					break;
				case 4:
					bindingEnergy = 0.95;
					break;
				case 5:
					bindingEnergy = 0.90;
					break;
				case 6:
					bindingEnergy = 0.90;
					break;
				case 7:
					bindingEnergy = 0.90;
					break;
				case 8:
					bindingEnergy = 0.85;
					break;
				case 9:
					bindingEnergy = 0.65;
					break;
				case 10:
					bindingEnergy = 0.65;
					break;
				case 11:
					bindingEnergy = 0.60;
					break;
				default:
					break;
				}
				break;
			case 8:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.05;
					break;
				case 2:
					bindingEnergy = 1.00;
					break;
				case 3:
					bindingEnergy = 0.95;
					break;
				case 4:
					bindingEnergy = 0.95;
					break;
				case 5:
					bindingEnergy = 0.90;
					break;
				case 6:
					bindingEnergy = 0.90;
					break;
				case 7:
					bindingEnergy = 0.85;
					break;
				case 8:
					bindingEnergy = 0.85;
					break;
				case 9:
					bindingEnergy = 0.65;
					break;
				case 10:
					bindingEnergy = 0.65;
					break;
				case 11:
					bindingEnergy = 0.60;
					break;
				default:
					break;
				}
				break;
			case 9:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.05;
					break;
				case 2:
					bindingEnergy = 1.00;
					break;
				case 3:
					bindingEnergy = 0.95;
					break;
				case 4:
					bindingEnergy = 0.95;
					break;
				case 5:
					bindingEnergy = 0.85;
					break;
				case 6:
					bindingEnergy = 0.85;
					break;
				case 7:
					bindingEnergy = 0.85;
					break;
				case 8:
					bindingEnergy = 0.85;
					break;
				case 9:
					bindingEnergy = 0.65;
					break;
				case 10:
					bindingEnergy = 0.65;
					break;
				case 11:
					bindingEnergy = 0.60;
					break;
				default:
					break;
				}
				break;
			case 10:
				switch (hSize) {
				case 1:
					bindingEnergy = 1.00;
					break;
				case 2:
					bindingEnergy = 0.95;
					break;
				case 3:
					bindingEnergy = 0.90;
					break;
				case 4:
					bindingEnergy = 0.90;
					break;
				case 5:
					bindingEnergy = 0.85;
					break;
				case 6:
					bindingEnergy = 0.85;
					break;
				case 7:
					bindingEnergy = 0.85;
					break;
				case 8:
					bindingEnergy = 0.80;
					break;
				case 9:
					bindingEnergy = 0.65;
					break;
				case 10:
					bindingEnergy = 0.60;
					break;
				case 11:
					bindingEnergy = 0.60;
					break;
				default:
					break;
				}
				break;
			case 11:
				switch (hSize) {
				case 1:
					bindingEnergy = 0.95;
					break;
				case 2:
					bindingEnergy = 0.95;
					break;
				case 3:
					bindingEnergy = 0.90;
					break;
				case 4:
					bindingEnergy = 0.90;
					break;
				case 5:
					bindingEnergy = 0.85;
					break;
				case 6:
					bindingEnergy = 0.85;
					break;
				case 7:
					bindingEnergy = 0.85;
					break;
				case 8:
					bindingEnergy = 0.80;
					break;
				case 9:
					bindingEnergy = 0.65;
					break;
				case 10:
					bindingEnergy = 0.60;
					break;
				case 11:
					bindingEnergy = 0.60;
					break;
				default:
					break;
				}
				break;
			case 12:
				switch (hSize) {
				case 1:
					bindingEnergy = 0.95;
					break;
				case 2:
					bindingEnergy = 0.90;
					break;
				case 3:
					bindingEnergy = 0.90;
					break;
				case 4:
					bindingEnergy = 0.85;
					break;
				case 5:
					bindingEnergy = 0.85;
					break;
				case 6:
					bindingEnergy = 0.85;
					break;
				case 7:
					bindingEnergy = 0.80;
					break;
				case 8:
					bindingEnergy = 0.80;
					break;
				case 9:
					bindingEnergy = 0.60;
					break;
				case 10:
					bindingEnergy = 0.60;
					break;
				case 11:
					bindingEnergy = 0.55;
					break;
				default:
					break;
				}
				break;
			case 13:
				switch (hSize) {
				case 1:
					bindingEnergy = 0.90;
					break;
				case 2:
					bindingEnergy = 0.90;
					break;
				case 3:
					bindingEnergy = 0.85;
					break;
				case 4:
					bindingEnergy = 0.85;
					break;
				case 5:
					bindingEnergy = 0.85;
					break;
				case 6:
					bindingEnergy = 0.85;
					break;
				case 7:
					bindingEnergy = 0.80;
					break;
				case 8:
					bindingEnergy = 0.80;
					break;
				case 9:
					bindingEnergy = 0.60;
					break;
				case 10:
					bindingEnergy = 0.60;
					break;
				case 11:
					bindingEnergy = 0.55;
					break;
				default:
					break;
				}
				break;
			case 14:
				switch (hSize) {
				case 1:
					bindingEnergy = 0.90;
					break;
				case 2:
					bindingEnergy = 0.90;
					break;
				case 3:
					bindingEnergy = 0.85;
					break;
				case 4:
					bindingEnergy = 0.85;
					break;
				case 5:
					bindingEnergy = 0.80;
					break;
				case 6:
					bindingEnergy = 0.80;
					break;
				case 7:
					bindingEnergy = 0.80;
					break;
				case 8:
					bindingEnergy = 0.70;
					break;
				case 9:
					bindingEnergy = 0.60;
					break;
				case 10:
					bindingEnergy = 0.60;
					break;
				case 11:
					bindingEnergy = 0.55;
					break;
				default:
					break;
				}
				break;
			}
		}
	}

//		std::cout << reaction.first.getName() << " + "
//				<< reaction.second.getName() << " <- "
//				<< reaction.dissociating.getName() << " : "
//				<< max(bindingEnergy, -5.0) << " = " << reaction.first.getFormationEnergy() << " + "
//				<< reaction.second.getFormationEnergy() << " - "
//				<< reaction.dissociating.getFormationEnergy() << std::endl;

	return max(bindingEnergy, -5.0);
}

IReactant * PSIClusterReactionNetwork::getSuperFromComp(IReactant::SizeType nHe,
		IReactant::SizeType nD, IReactant::SizeType nT,
		IReactant::SizeType nV) {

	// Requests for finding a particular supercluster have high locality.
	// See if the last supercluster we were asked to find is the right
	// one for this request.
	static IReactant* lastRet;
	if (lastRet
			and static_cast<PSISuperCluster*>(lastRet)->isIn(nHe, nD, nT, nV)) {
		return lastRet;
	}

	// We didn't find the last supercluster in our cache, so do a full lookup.
	IReactant* ret = nullptr;

	for (auto const& superMapItem : getAll(ReactantType::PSISuper)) {

		auto const& reactant =
				static_cast<PSISuperCluster&>(*(superMapItem.second));
		if (reactant.isIn(nHe, nD, nT, nV)) {
			lastRet = superMapItem.second.get();
			return superMapItem.second.get();
		}
	}

	return ret;
}

} // namespace xolotlCore
