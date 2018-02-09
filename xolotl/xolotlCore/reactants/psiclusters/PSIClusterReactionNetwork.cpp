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
				ReactantType::HeV, ReactantType::HeI, ReactantType::PSISuper },
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
	auto dissociationAllowed = canDissociate(prref);

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

void PSIClusterReactionNetwork::createReactionConnectivity() {
	// Initial declarations
	IReactant::SizeType firstSize = 0, secondSize = 0, productSize = 0;

	// Single species clustering (He, V, I)
	// We know here that only Xe_1 can cluster so we simplify the search
	// X_(a-i) + X_i --> X_a
	// Make a vector of types
	std::vector<ReactantType> typeVec { ReactantType::He, ReactantType::V,
			ReactantType::I };
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

	// Helium absorption by HeV clusters
	// He_(a) + (He_b)(V_c) --> [He_(a+b)](V_c)
	// Consider each He reactant.
	for (auto const& heMapItem : getAll(ReactantType::He)) {

		auto& heReactant = *(heMapItem.second);

		// Skip if it can't diffuse
		if (xolotlCore::equal(heReactant.getDiffusionFactor(), 0.0))
			continue;
		// Get its size
		firstSize = heReactant.getSize();

		// Consider product with each HeV cluster
		for (auto const& heVMapItem : getAll(ReactantType::HeV)) {

			auto& heVReactant = *(heVMapItem.second);

			// Get its composition
			auto& comp = heVReactant.getComposition();
			// Create the composition of the potential product
			auto newNumHe = comp[toCompIdx(Species::He)] + firstSize;
			auto newNumV = comp[toCompIdx(Species::V)];

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::HeV, newComp);

			// Check if the product can be a super cluster
			if (!product) {
				// Check if it is a super cluster from the map
				product = getSuperFromComp(newNumHe, newNumV);
			}
			// Check that the reaction can occur
			if (product
					&& (heReactant.getDiffusionFactor() > 0.0
							|| heVReactant.getDiffusionFactor() > 0.0)) {

				defineProductionReaction(heReactant, heVReactant, *product,
						newNumHe, newNumV);
			}
		}

		// Consider product with each super cluster
		for (auto const& superMapItem : getAll(ReactantType::PSISuper)) {

			auto& superCluster =
					static_cast<PSISuperCluster&>(*(superMapItem.second));
			std::vector<PendingProductionReactionInfo> prInfos;

			// Get its boundaries
			auto const& vBounds = superCluster.getVBounds();
			// Loop on them
			for (auto const& i : superCluster.getHeBounds()) {
				for (auto const& j : superCluster.getVBounds()) {
					// Assume the product can only be a super cluster here
					auto newNumHe = i + firstSize;
					auto newNumV = j;
					IReactant* product = getSuperFromComp(newNumHe, newNumV);
					// Check that the reaction can occur
					if (product
							&& (heReactant.getDiffusionFactor() > 0.0
									|| superCluster.getDiffusionFactor() > 0.0)) {

						// Note that current reactant reacts with
						// current superCluster to produce product,
						// according to current parameters.
						prInfos.emplace_back(*product, newNumHe, newNumV, i, j);
					}
				}
			}

			// Now that we know how current reactant reacts with
			// current superCluster, create the production
			// reaction(s) for them.
			defineProductionReactions(heReactant, superCluster, prInfos);
		}
	}

	// Vacancy absorption by HeV clusters
	// (He_a)(V_b) + V_c --> (He_a)[V_(b+c)]
	// Consider each V cluster.
	for (auto const& vMapItem : getAll(ReactantType::V)) {

		auto& vReactant = *(vMapItem.second);

		// Skip if it can't diffuse
		if (xolotlCore::equal(vReactant.getDiffusionFactor(), 0.0))
			continue;
		// Get the V size
		firstSize = vReactant.getSize();
		// Consider product with every HeV cluster.
		for (auto const& heVMapItem : getAll(ReactantType::HeV)) {

			auto& heVReactant = *(heVMapItem.second);

			// Get its composition
			auto& comp = heVReactant.getComposition();
			// Create the composition of the potential product
			auto newNumHe = comp[toCompIdx(Species::He)];
			auto newNumV = comp[toCompIdx(Species::V)] + firstSize;

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::HeV, newComp);

			// Check if the product can be a super cluster
			if (!product) {
				product = getSuperFromComp(newNumHe, newNumV);
			}
			// Check that the reaction can occur
			if (product
					&& (vReactant.getDiffusionFactor() > 0.0
							|| heVReactant.getDiffusionFactor() > 0.0)) {

				defineProductionReaction(vReactant, heVReactant, *product,
						newNumHe, newNumV);
			}
		}

		// Consider product with super clusters.
		for (auto const& superMapItem : getAll(ReactantType::PSISuper)) {

			auto& superCluster =
					static_cast<PSISuperCluster&>(*(superMapItem.second));
			std::vector<PendingProductionReactionInfo> prInfos;

			// Get its boundaries
			auto const& heBounds = superCluster.getHeBounds();
			auto const& vBounds = superCluster.getVBounds();
			// Loop on them
			for (auto const& i : superCluster.getHeBounds()) {
				for (auto const& j : superCluster.getVBounds()) {

					// Assume the product can only be a super cluster here
					auto newNumHe = i;
					auto newNumV = j + firstSize;
					IReactant* product = getSuperFromComp(newNumHe, newNumV);
					// Check that the reaction can occur
					if (product
							&& (vReactant.getDiffusionFactor() > 0.0
									|| superCluster.getDiffusionFactor() > 0.0)) {
						prInfos.emplace_back(*product, newNumHe, newNumV, i, j);
					}
				}
			}

			// Now that we know how current reactant interacts with
			// current supercluster, define the production reactions.
			defineProductionReactions(vReactant, superCluster, prInfos);
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
			auto newNumHe = firstSize;
			auto newNumV = secondSize;

			// Get the product
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::HeV, newComp);

			// Check if the product can be a super cluster
			if (!product) {
				product = getSuperFromComp(newNumHe, newNumV);
			}
			// Check that the reaction can occur
			if (product
					&& (heReactant.getDiffusionFactor() > 0.0
							|| vReactant.getDiffusionFactor() > 0.0)) {

				defineProductionReaction(heReactant, vReactant, *product,
						newNumHe, newNumV);
			}
		}
	}

	// Vacancy reduction by Interstitial absorption in HeV clusters
	// (He_a)(V_b) + (I_c) --> (He_a)[V_(b-c)]
	// Consider each I cluster
	for (auto const& iMapItem : getAll(ReactantType::I)) {

		auto& iReactant = *(iMapItem.second);

		// Get its size
		firstSize = iReactant.getSize();
		// Consider product with each HeV cluster.
		for (auto const& heVMapItem : getAll(ReactantType::HeV)) {

			auto& heVReactant = *(heVMapItem.second);

			// Get its composition
			auto& comp = heVReactant.getComposition();
			// The product can be He or HeV
			IReactant * product = nullptr;
			if (comp[toCompIdx(Species::V)] == firstSize) {
				// The product is He
				product = get(Species::He, comp[toCompIdx(Species::He)]);
			} else {
				// The product is HeV
				// Create the composition of the potential product
				IReactant::Composition newComp;
				newComp[toCompIdx(Species::He)] = comp[toCompIdx(Species::He)];
				newComp[toCompIdx(Species::V)] = comp[toCompIdx(Species::V)]
						- firstSize;
				// Get the product
				product = get(ReactantType::HeV, newComp);
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

			// Get its boundaries
			auto const& heBounds = superCluster.getHeBounds();
			auto const& vBounds = superCluster.getVBounds();
			// Loop on them
			for (auto const& i : superCluster.getHeBounds()) {
				for (auto const& j : superCluster.getVBounds()) {
					// The product might be HeV or He
					auto newNumHe = i;
					auto newNumV = j - firstSize;

					// Get the product
					IReactant* product = nullptr;
					if (newNumV == 0) {
						// The product is He
						product = get(Species::He, i);
					} else {
						// Create the composition of the potential product
						IReactant::Composition newComp;
						newComp[toCompIdx(Species::He)] = newNumHe;
						newComp[toCompIdx(Species::V)] = newNumV;
						product = get(ReactantType::HeV, newComp);

						// If the product doesn't exist check for super clusters
						if (!product) {
							product = getSuperFromComp(newNumHe, newNumV);
						}
					}
					// Check that the reaction can occur
					if (product
							&& (iReactant.getDiffusionFactor() > 0.0
									|| superCluster.getDiffusionFactor() > 0.0)) {

						prInfos.emplace_back(*product, newNumHe, newNumV, i, j);
					}
				}
			}

			// Now that we know how current reactant interacts with
			// current supercluster, define its production reactions
			// according to given parameters.
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
			for (auto const& iMapItem : getAll(ReactantType::I)) {
				auto& iReactant = *(iMapItem.second);
				// Get the size of the I cluster
				int iSize = iReactant.getSize();
				// Create the composition of the potential product
				IReactant::Composition newComp;
				int newNumHe = firstSize + secondSize;
				int newNumV = iSize;
				newComp[toCompIdx(Species::He)] = newNumHe;
				newComp[toCompIdx(Species::V)] = newNumV;
				product = get(ReactantType::HeV, newComp);

				// Check if the product can be a super cluster
				if (!product) {
					// Check if it is a super cluster from the map
					product = getSuperFromComp(newNumHe, newNumV);
				}

				// Check that the reaction can occur
				if (product
						&& (firstReactant.getDiffusionFactor() > 0.0
								|| secondReactant.getDiffusionFactor() > 0.0)) {

					defineProductionReaction(firstReactant, secondReactant,
							*product, newNumHe, newNumV);
					// This is a reaction with two products so we need to tell the other product
					// it is participating too
					defineProductionReaction(firstReactant, secondReactant,
							iReactant, newNumHe, newNumV, 0, 0, true);

					// Stop the loop on I clusters here
					break;
				}
			}
		}
	}

	// Helium absorption by HeV leading to trap mutation
	// (He_a)(V_b) + He_c --> [He_(a+c)][V_(b+d)] + I_d
	// Loop on the He clusters
	for (auto const& heMapItem : getAll(ReactantType::He)) {
		auto& heReactant = *(heMapItem.second);

		// Skip if it can't diffuse
		if (xolotlCore::equal(heReactant.getDiffusionFactor(), 0.0))
			continue;
		// Get its size
		firstSize = heReactant.getSize();

		// Loop on the HeV clusters
		for (auto const& heVMapItem : getAll(ReactantType::HeV)) {

			auto& heVReactant = *(heVMapItem.second);

			// Get its composition
			auto& comp = heVReactant.getComposition();
			// Create the composition of the potential product
			auto newNumHe = comp[toCompIdx(Species::He)] + firstSize;
			auto newNumV = comp[toCompIdx(Species::V)];

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::V)] = newNumV;
			auto product = get(ReactantType::HeV, newComp);
			// Doesn't do anything if the product exist
			if (product)
				continue;

			// Trap mutation is happening
			// Loop on the possible I starting by the smallest
			for (auto const& iMapItem : getAll(ReactantType::I)) {
				auto& iReactant = *(iMapItem.second);
				// Get the size of the I cluster
				int iSize = iReactant.getSize();
				// Create the composition of the potential product
				newNumV += iSize;
				newComp[toCompIdx(Species::V)] = newNumV;
				product = get(ReactantType::HeV, newComp);

				// Check if the product can be a super cluster
				if (!product) {
					// Check if it is a super cluster from the map
					product = getSuperFromComp(newNumHe, newNumV);
				}

				// Check that the reaction can occur
				if (product
						&& (heReactant.getDiffusionFactor() > 0.0
								|| heVReactant.getDiffusionFactor() > 0.0)) {

					defineProductionReaction(heReactant, heVReactant, *product,
							newNumHe, newNumV);
					// This is a reaction with two products so we need to tell the other product
					// it is participating too
					// Define the basic production reaction.
					auto& reaction = defineReactionBase(heReactant,
							heVReactant);
					iReactant.resultFrom(reaction, newNumHe, newNumV);

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

			// Get its boundaries
			auto const& heBounds = superCluster.getHeBounds();
			auto const& vBounds = superCluster.getVBounds();
			// Loop on them
			for (auto const& i : superCluster.getHeBounds()) {
				for (auto const& j : superCluster.getVBounds()) {
					// Check these coordinates are actually contained by the super cluster
					if (!superCluster.isIn(i, j))
						continue;

					// The product might be HeV or He
					auto newNumHe = i + firstSize;
					auto newNumV = j;

					// Get the product
					IReactant* product = getSuperFromComp(newNumHe, newNumV);
					// Skip if the product exists because we want trap mutation
					if (product)
						continue;

					// Trap mutation is happening
					// Loop on the possible I starting by the smallest
					for (auto const& iMapItem : getAll(ReactantType::I)) {
						auto& iReactant = *(iMapItem.second);
						// Get the size of the I cluster
						int iSize = iReactant.getSize();
						// Update the composition of the potential product
						newNumV += iSize;
						product = getSuperFromComp(newNumHe, newNumV);

						// Check that the reaction can occur
						if (product && heReactant.getDiffusionFactor() > 0.0) {

							prInfos1.emplace_back(*product, newNumHe, newNumV,
									i, j);

							prInfos2.emplace_back(iReactant, 0, 0, i, j);
						}
					}
				}
			}

			// Now that we know how current reactant interacts with
			// current supercluster, define its production reactions
			// according to given parameters.
			defineProductionReactions(heReactant, superCluster, prInfos1);
			defineProductionReactions(heReactant, superCluster, prInfos2, true);
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
			auto newNumHe = comp[toCompIdx(Species::He)] + firstSize;
			auto newNumI = comp[toCompIdx(Species::I)];

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::I)] = newNumI;
			auto product = get(ReactantType::HeI, newComp);

			// Check that the reaction can occur
			if (product
					&& (heReactant.getDiffusionFactor() > 0.0
							|| heIReactant.getDiffusionFactor() > 0.0)) {

				defineProductionReaction(heReactant, heIReactant, *product,
						newNumHe, newNumI);
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
			auto newNumHe = comp[toCompIdx(Species::He)];
			auto newNumI = comp[toCompIdx(Species::I)] + firstSize;

			// Check if product already exists.
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::I)] = newNumI;
			auto product = get(ReactantType::HeI, newComp);

			// Check that the reaction can occur
			if (product
					&& (iReactant.getDiffusionFactor() > 0.0
							|| heIReactant.getDiffusionFactor() > 0.0)) {

				defineProductionReaction(iReactant, heIReactant, *product,
						newNumHe, newNumI);
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
			auto newNumHe = firstSize;
			auto newNumI = secondSize;

			// Get the product
			IReactant::Composition newComp;
			newComp[toCompIdx(Species::He)] = newNumHe;
			newComp[toCompIdx(Species::I)] = newNumI;
			auto product = get(ReactantType::HeI, newComp);

			// Check that the reaction can occur
			if (product
					&& (heReactant.getDiffusionFactor() > 0.0
							|| iReactant.getDiffusionFactor() > 0.0)) {

				defineProductionReaction(heReactant, iReactant, *product,
						newNumHe, newNumI);
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

bool PSIClusterReactionNetwork::canDissociate(
		ProductionReaction& reaction) const {
	// Assume reaction can dissociate by default.
	bool ret = true;

	// Check if at least one of the potentially emitted cluster is size one
	if (reaction.first.getSize() != 1 && reaction.second.getSize() != 1) {
		// Don't add the reverse reaction
		ret = false;
	}

	// Check for trap mutations (with XOR)
	else if ((reaction.first.getType() == ReactantType::I)
			== !(reaction.second.getType() == ReactantType::I)) {
		// Don't add the reverse reaction
		ret = false;
	}

	return ret;
}

void PSIClusterReactionNetwork::checkForDissociation(
		IReactant& emittingReactant, ProductionReaction& reaction, int a, int b,
		int c, int d) {

	// Check if reaction can dissociate.
	if (canDissociate(reaction)) {
		// The dissociation can occur, so create a reaction for it.
		defineDissociationReaction(reaction, emittingReactant, a, b, c, d);
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
				currReactant.setHeMomentumId(id);
				currReactant.setVMomentumId(id);
			});

	// Get all the super clusters and loop on them
	for (auto const& currMapItem : clusterTypeMap[ReactantType::PSISuper]) {

		auto& currCluster = static_cast<PSISuperCluster&>(*(currMapItem.second));

		id++;
		currCluster.setHeMomentumId(id);
		id++;
		currCluster.setVMomentumId(id);

		// Update the HeV size
		IReactant::SizeType clusterSize = currCluster.getHeBounds().second
				+ currCluster.getVBounds().second;
		if (clusterSize > maxClusterSizeMap[ReactantType::HeV]) {
			maxClusterSizeMap[ReactantType::HeV] = clusterSize;
		}
	}

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
			[&concentrations](const ReactantMap::value_type& currMapItem) {

				auto& cluster = static_cast<PSISuperCluster&>(*(currMapItem.second));

				cluster.setZerothMomentum(concentrations[cluster.getId() - 1]);
				cluster.setHeMomentum(concentrations[cluster.getHeMomentumId() - 1]);
				cluster.setVMomentum(concentrations[cluster.getVMomentumId() - 1]);
			});

	return;
}

void PSIClusterReactionNetwork::getDiagonalFill(int *diagFill) {
	// Degrees of freedom is the total number of clusters in the network
	const int dof = getDOF();

	// Get the connectivity for each reactant
	std::for_each(allReactants.begin(), allReactants.end(),
			[&diagFill,&dof,this](const IReactant& reactant) {

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
	for (auto const& superMapItem : getAll(ReactantType::PSISuper)) {

		auto const& reactant =
				static_cast<PSISuperCluster&>(*(superMapItem.second));

		// Get the reactant and its connectivity
		auto const& connectivity = reactant.getConnectivity();
		auto connectivityLength = connectivity.size();
		// Get the helium momentum id so that the connectivity can be lined up in
		// the proper column
		auto id = reactant.getHeMomentumId() - 1;

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

		// Get the vacancy momentum id so that the connectivity can be lined up in
		// the proper column
		id = reactant.getVMomentumId() - 1;

		// Add it to the diagonal fill block
		for (int j = 0; j < connectivityLength; j++) {
			// The id starts at j*connectivity length and is always offset
			// by the id, which denotes the exact column.
			auto index = (id) * dof + j;
			diagFill[index] = connectivity[j];
		}
		// Update the map
		dFillMap[id] = columnIds;
	}

	return;
}

double PSIClusterReactionNetwork::getTotalAtomConcentration() {
	// Initial declarations
	double heliumConc = 0.0;

	// Sum over all He clusters.
	for (auto const& currMapItem : getAll(ReactantType::He)) {

		// Get the cluster and its composition
		auto const& cluster = *(currMapItem.second);
		double size = cluster.getSize();

		// Add the concentration times the He content to the total helium concentration
		heliumConc += cluster.getConcentration() * size;
	}

	// Sum over all HeV clusters.
	for (auto const& currMapItem : getAll(ReactantType::HeV)) {

		// Get the cluster and its composition
		auto const& cluster = *(currMapItem.second);
		auto& comp = cluster.getComposition();

		// Add the concentration times the He content to the total helium concentration
		heliumConc += cluster.getConcentration() * comp[toCompIdx(Species::He)];
	}

	// Sum over all super clusters.
	for (auto const& currMapItem : getAll(ReactantType::PSISuper)) {

		// Get the cluster
		auto const& cluster =
				static_cast<PSISuperCluster&>(*(currMapItem.second));

		// Add its total helium concentration helium concentration
		heliumConc += cluster.getTotalHeliumConcentration();
	}

	return heliumConc;
}

double PSIClusterReactionNetwork::getTotalTrappedAtomConcentration() {
	// Initial declarations
	double heliumConc = 0.0;

	// Sum over all HeV clusters.
	for (auto const& currMapItem : getAll(ReactantType::HeV)) {
		// Get the cluster and its composition
		auto const& cluster = *(currMapItem.second);
		auto& comp = cluster.getComposition();

		// Add the concentration times the He content to the total helium concentration
		heliumConc += cluster.getConcentration() * comp[toCompIdx(Species::He)];
	}

	// Sum over all super clusters.
	for (auto const& currMapItem : getAll(ReactantType::PSISuper)) {
		// Get the cluster
		auto const& cluster =
				static_cast<PSISuperCluster&>(*(currMapItem.second));

		// Add its total helium concentration
		heliumConc += cluster.getTotalHeliumConcentration();
	}

	return heliumConc;
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
	for (auto const& currMapItem : getAll(ReactantType::HeV)) {
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

		// Compute the helium momentum flux
		auto flux = superCluster.getHeMomentumFlux();
		// Update the concentration of the cluster
		auto reactantIndex = superCluster.getHeMomentumId() - 1;
		updatedConcOffset[reactantIndex] += flux;

		// Compute the vacancy momentum flux
		flux = superCluster.getVMomentumFlux();
		// Update the concentration of the cluster
		reactantIndex = superCluster.getVMomentumId() - 1;
		updatedConcOffset[reactantIndex] += flux;
	}

	return;
}

void PSIClusterReactionNetwork::computeAllPartials(double *vals, int *indices,
		int *size) {
	// Initial declarations
	const int dof = getDOF();
	std::vector<double> clusterPartials;
	clusterPartials.resize(dof, 0.0);

	// Get the super clusters
	auto const& superClusters = getAll(ReactantType::PSISuper);

	// Make a vector of types for the non super clusters
	std::vector<ReactantType> typeVec { ReactantType::He, ReactantType::V,
			ReactantType::I, ReactantType::HeV };
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
			// Number of partial derivatives
			auto pdColIdsVectorSize = pdColIdsVector.size();
			size[reactantIndex] = pdColIdsVectorSize;

			// Loop over the list of column ids
			for (int j = 0; j < pdColIdsVectorSize; j++) {
				// Set the index
				indices[reactantIndex * dof + j] = pdColIdsVector[j];

				// Get the partial derivative from the array of all of the partials
				vals[reactantIndex * dof + j] =
						clusterPartials[pdColIdsVector[j]];

				// Reset the cluster partial value to zero. This is much faster
				// than using memset.
				clusterPartials[pdColIdsVector[j]] = 0.0;
			}
		}
	}

	// Update the column in the Jacobian that represents the moment for the super clusters
	for (auto const& currMapItem : superClusters) {

		auto const& reactant =
				static_cast<PSISuperCluster&>(*(currMapItem.second));

		{
			// Get the super cluster index
			auto reactantIndex = reactant.getId() - 1;

			// Get the partial derivatives
			reactant.getPartialDerivatives(clusterPartials);

			// Get the list of column ids from the map
			auto const& pdColIdsVector = dFillMap.at(reactantIndex);
			// Number of partial derivatives
			auto pdColIdsVectorSize = pdColIdsVector.size();
			size[reactantIndex] = pdColIdsVectorSize;

			// Loop over the list of column ids
			for (int j = 0; j < pdColIdsVectorSize; j++) {
				// Set the index
				indices[reactantIndex * dof + j] = pdColIdsVector[j];
				// Get the partial derivative from the array of all of the partials
				vals[reactantIndex * dof + j] =
						clusterPartials[pdColIdsVector[j]];

				// Reset the cluster partial value to zero. This is much faster
				// than using memset.
				clusterPartials[pdColIdsVector[j]] = 0.0;
			}
		}

		{
			// Get the helium momentum index
			auto reactantIndex = reactant.getHeMomentumId() - 1;

			// Get the partial derivatives
			reactant.getHeMomentPartialDerivatives(clusterPartials);
			// Get the list of column ids from the map
			auto const& pdColIdsVector = dFillMap.at(reactantIndex);
			// Number of partial derivatives
			auto pdColIdsVectorSize = pdColIdsVector.size();
			size[reactantIndex] = pdColIdsVectorSize;

			// Loop over the list of column ids
			for (int j = 0; j < pdColIdsVectorSize; j++) {
				// Set the index
				indices[reactantIndex * dof + j] = pdColIdsVector[j];
				// Get the partial derivative from the array of all of the partials
				vals[reactantIndex * dof + j] =
						clusterPartials[pdColIdsVector[j]];

				// Reset the cluster partial value to zero. This is much faster
				// than using memset.
				clusterPartials[pdColIdsVector[j]] = 0.0;
			}
		}

		{
			// Get the vacancy momentum index
			auto reactantIndex = reactant.getVMomentumId() - 1;

			// Get the partial derivatives
			reactant.getVMomentPartialDerivatives(clusterPartials);
			// Get the list of column ids from the map
			auto const& pdColIdsVector = dFillMap.at(reactantIndex);
			// Number of partial derivatives
			auto pdColIdsVectorSize = pdColIdsVector.size();
			size[reactantIndex] = pdColIdsVectorSize;

			// Loop over the list of column ids
			for (int j = 0; j < pdColIdsVectorSize; j++) {
				// Set the index
				indices[reactantIndex * dof + j] = pdColIdsVector[j];
				// Get the partial derivative from the array of all of the partials
				vals[reactantIndex * dof + j] =
						clusterPartials[pdColIdsVector[j]];

				// Reset the cluster partial value to zero. This is much faster
				// than using memset.
				clusterPartials[pdColIdsVector[j]] = 0.0;
			}
		}
	}

	return;
}

double PSIClusterReactionNetwork::computeBindingEnergy(
		const DissociationReaction& reaction) const {
	// for the dissociation A --> B + C we need A binding energy
	// E_b(A) = E_f(B) + E_f(C) - E_f(A) where E_f is the formation energy
	double bindingEnergy = reaction.first.getFormationEnergy()
			+ reaction.second.getFormationEnergy()
			- reaction.dissociating.getFormationEnergy();

	return max(bindingEnergy, -5.0);
}

IReactant * PSIClusterReactionNetwork::getSuperFromComp(IReactant::SizeType nHe,
		IReactant::SizeType nV) {

	// Requests for finding a particular supercluster have high locality.
	// See if the last supercluster we were asked to find is the right
	// one for this request.
	static IReactant* lastRet = nullptr;
	if (lastRet and static_cast<PSISuperCluster*>(lastRet)->isIn(nHe, nV)) {
		return lastRet;
	}

	// We didn't find the last supercluster in our cache, so do a full lookup.
	IReactant* ret = nullptr;

	for (auto const& superMapItem : getAll(ReactantType::PSISuper)) {

		auto const& reactant =
				static_cast<PSISuperCluster&>(*(superMapItem.second));
		if (reactant.isIn(nHe, nV))
			return superMapItem.second.get();
	}

	return ret;
}

} // namespace xolotlCore

