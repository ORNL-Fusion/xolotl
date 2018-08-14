#include <cassert>
#include <iterator>
#include "FeClusterReactionNetwork.h"
#include "FeCluster.h"
#include "FeSuperCluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

namespace xolotlCore {

FeClusterReactionNetwork::FeClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork( { ReactantType::V, ReactantType::I, ReactantType::He,
				ReactantType::HeV, ReactantType::FeSuper },
				ReactantType::FeSuper, registry) {

	// Initialize default properties
	dissociationsEnabled = true;

	return;
}

double FeClusterReactionNetwork::calculateDissociationConstant(
		const DissociationReaction& reaction, int i) const {

	// If the dissociations are not allowed
	if (!dissociationsEnabled)
		return 0.0;

	// The atomic volume is computed by considering the BCC structure of the
	// iron. In a given lattice cell in iron there are iron atoms
	// at each corner and a iron atom in the center. The iron atoms at
	// the corners are shared across a total of eight cells. The fraction of
	// the volume of the lattice cell that is filled with iron atoms is the
	// atomic volume and is a_0^3/(8*1/8 + 1) = 0.5*a_0^3.
	double atomicVolume = 0.5 * xolotlCore::ironLatticeConstant
			* xolotlCore::ironLatticeConstant * xolotlCore::ironLatticeConstant;

	// Get the rate constant from the reverse reaction
	double kPlus = reaction.reverseReaction->kConstant[i];

	// Calculate and return
	double bindingEnergy = computeBindingEnergy(reaction);
	double k_minus_exp = exp(
			-1.0 * bindingEnergy / (xolotlCore::kBoltzmann * temperature));
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;

	return k_minus;
}

void FeClusterReactionNetwork::defineProductionReactions(IReactant& r1,
		IReactant& r2, IReactant& product) {

	// Check if the reaction can happen
	auto const& superProd = static_cast<FeCluster const&>(product);
	auto const& heBounds = superProd.getHeBounds();
	auto const& vBounds = superProd.getVBounds();
	int productLoHe = *(heBounds.begin()), productHiHe = *(heBounds.end()) - 1,
			productLoV = *(vBounds.begin()), productHiV = *(vBounds.end()) - 1,
			loHe = 0, hiHe = 0, loV = 0, hiV = 0, singleHeSize = 0,
			singleVSize = 0;

	if (r1.getType() == ReactantType::FeSuper) {
		auto const& super = static_cast<FeSuperCluster const&>(r1);
		auto const& heBounds = super.getHeBounds();
		loHe = *(heBounds.begin());
		hiHe = *(heBounds.end()) - 1;
		auto const& vBounds = super.getVBounds();
		loV = *(vBounds.begin());
		hiV = *(vBounds.end()) - 1;
		auto singleComp = r2.getComposition();
		singleHeSize = singleComp[toCompIdx(Species::He)];
		singleVSize = singleComp[toCompIdx(Species::V)]
				- singleComp[toCompIdx(Species::I)]; // can be < 0
	}
	if (r2.getType() == ReactantType::FeSuper) {
		auto const& super = static_cast<FeSuperCluster const&>(r2);
		auto const& heBounds = super.getHeBounds();
		loHe = *(heBounds.begin());
		hiHe = *(heBounds.end()) - 1;
		auto const& vBounds = super.getVBounds();
		loV = *(vBounds.begin());
		hiV = *(vBounds.end()) - 1;
		auto singleComp = r1.getComposition();
		singleHeSize = singleComp[toCompIdx(Species::He)];
		singleVSize = singleComp[toCompIdx(Species::V)]
				- singleComp[toCompIdx(Species::I)]; // can be < 0
	}

	int heWidth = std::min(productHiHe, hiHe + singleHeSize)
			- std::max(productLoHe, loHe + singleHeSize) + 1;
	int vWidth = std::min(productHiV, hiV + singleVSize)
			- std::max(productLoV, loV + singleVSize) + 1;

	if (heWidth <= 0)
		return;
	if (vWidth <= 0)
		return;

	// Define the production reaction to the network.
	// Do this *once* for the given reactants, since it doesn't
	// depend on the product or the other parameters.
	std::unique_ptr<ProductionReaction> reaction(
			new ProductionReaction(r1, r2));
	auto& prref = add(std::move(reaction));

	// Determine if reverse reaction is allowed.
	auto dissociationAllowed = canDissociate(prref);

	// Tell both reactants they are involved in this reaction
	// with the given parameters.
	prref.first.participateIn(prref, product);
	prref.second.participateIn(prref, product);
	product.resultFrom(prref, product);

	// Determine if reverse reaction is allowed.
	if (dissociationAllowed) {
		// Forward reaction can dissociate.
		// Define all dissociation reactions for this forward reaction
		defineDissociationReactions(prref, product);
	}

	return;
}

void FeClusterReactionNetwork::defineDissociationReactions(
		ProductionReaction& forwardReaction, IReactant& disso) {
	// Add a dissociation reaction to our network.
	// Do this once here for each forward reaction product.
	std::unique_ptr<DissociationReaction> dissociationReaction(
			new DissociationReaction(disso, forwardReaction.first,
					forwardReaction.second, &forwardReaction));
	auto& drref = add(std::move(dissociationReaction));

	// Tell all participants in this reaction of their involvement.
	drref.first.participateIn(drref, disso);
	drref.second.participateIn(drref, disso);
	disso.emitFrom(drref, disso);

	return;
}

void FeClusterReactionNetwork::createReactionConnectivity() {
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
			int newNumHe = comp[toCompIdx(Species::He)] + firstSize;
			int newNumV = comp[toCompIdx(Species::V)];

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

				int a[4] = { newNumHe, newNumV, 0, 0 };
				defineProductionReaction(heReactant, heVReactant, *product, a);
			}
		}

		// Consider product with each super cluster
		for (auto const& superMapItem : getAll(ReactantType::FeSuper)) {
			// Loop on the potential super products
			for (auto const& superMapItemProd : getAll(ReactantType::FeSuper)) {

				// This method will check if the reaction is possible and then add it to the list
				defineProductionReactions(heReactant, *(superMapItem.second),
						*(superMapItemProd.second));
			}
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
			int newNumHe = comp[toCompIdx(Species::He)];
			int newNumV = comp[toCompIdx(Species::V)] + firstSize;

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

				int a[4] = { newNumHe, newNumV, 0, 0 };
				defineProductionReaction(vReactant, heVReactant, *product, a);
			}
		}

		// Consider product with each super cluster
		for (auto const& superMapItem : getAll(ReactantType::FeSuper)) {
			// Loop on the potential super products
			for (auto const& superMapItemProd : getAll(ReactantType::FeSuper)) {

				// This method will check if the reaction is possible and then add it to the list
				defineProductionReactions(vReactant, *(superMapItem.second),
						*(superMapItemProd.second));
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
			auto product = get(ReactantType::HeV, newComp);

			// Check if the product can be a super cluster
			if (!product) {
				product = getSuperFromComp(newNumHe, newNumV);
			}
			// Check that the reaction can occur
			if (product
					&& (heReactant.getDiffusionFactor() > 0.0
							|| vReactant.getDiffusionFactor() > 0.0)) {

				int a[4] = { newNumHe, newNumV, 0, 0 };
				defineProductionReaction(heReactant, vReactant, *product, a);
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

		// Consider product with each super cluster
		for (auto const& superMapItem : getAll(ReactantType::FeSuper)) {
			// Loop on the potential super products
			for (auto const& superMapItemProd : getAll(ReactantType::FeSuper)) {

				// This method will check if the reaction is possible and then add it to the list
				defineProductionReactions(iReactant, *(superMapItem.second),
						*(superMapItemProd.second));
			}
			// Loop on the potential He products
			for (auto const& heMapItemProd : getAll(ReactantType::He)) {

				// This method will check if the reaction is possible and then add it to the list
				defineProductionReactions(iReactant, *(superMapItem.second),
						*(heMapItemProd.second));
			}
			// Loop on the potential HeV products
			for (auto const& heVMapItemProd : getAll(ReactantType::HeV)) {

				// This method will check if the reaction is possible and then add it to the list
				defineProductionReactions(iReactant, *(superMapItem.second),
						*(heVMapItemProd.second));
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

	return;
}

bool FeClusterReactionNetwork::canDissociate(
		ProductionReaction& reaction) const {
	// Assume reaction can dissociate by default.
	bool ret = true;

	// Check if at least one of the potentially emitted cluster is size one
	if (reaction.first.getSize() != 1 && reaction.second.getSize() != 1) {
		// Don't add the reverse reaction
		ret = false;
	}
	// remove He+He
	else if (reaction.first.getSize() == 1 && reaction.second.getSize() == 1
			&& reaction.first.getType() == ReactantType::He
			&& reaction.second.getType() == ReactantType::He) {
		// Don't add the reverse reaction
		ret = false;
	}

//	// Check for trap mutations (with XOR)
//	else if ((reaction->first.getType() == ReactantType::I)
//			== !(reaction->second.getType() == ReactantType::I)) {
//		// Don't add the reverse reaction
//		return;
//	}

	return ret;
}

void FeClusterReactionNetwork::checkForDissociation(IReactant& emittingReactant,
		ProductionReaction& reaction, int a[4], int b[4]) {

	// Check if reaction can dissociate.
	if (canDissociate(reaction)) {
		// The dissociation can occur, so create a reaction for it.
		defineDissociationReaction(reaction, emittingReactant, a, b);
	}

	return;
}

void FeClusterReactionNetwork::setTemperature(double temp, int i) {
	ReactionNetwork::setTemperature(temp, i);

	computeRateConstants(i);

	return;
}

void FeClusterReactionNetwork::buildSuperClusterMap(
		const std::vector<IReactant::SizeType>& bounds) {

	// Save the bounds to use.
	boundVector = bounds;

	// Since we represent the super cluster lookup map using a dense
	// data structure, we must initialize every entry with something
	// that signifies 'Invalid.'
	auto const& superClusters = clusterTypeMap[ReactantType::FeSuper];
	auto bvSize = boundVector.size();
	superClusterLookupMap.resize(bvSize);
	std::for_each(superClusterLookupMap.begin(), superClusterLookupMap.end(),
			[bvSize,&superClusters](HeVToSuperClusterMap::value_type& currVector) {

				currVector.resize(bvSize, superClusters.end());
			});

	// Build a map of super clusters, keyed by (baseHe, baseV) pairs
	// where base* indicates the lower bound of the super cluster's
	// interval for that species type.
	for (auto iter = superClusters.begin(); iter != superClusters.end();
			++iter) {
		// Add the super cluster to our lookup map based on
		// its He and V intervals.
		auto& currCluster = *(iter->second);
		auto const& comp = currCluster.getComposition();
		auto currHe = comp[toCompIdx(Species::He)];
		auto currV = comp[toCompIdx(Species::V)];
		auto heIntervalIdx = findBoundsIntervalBaseIdx(currHe);
		auto vIntervalIdx = findBoundsIntervalBaseIdx(currV);

		superClusterLookupMap[heIntervalIdx][vIntervalIdx] = iter;
	}
}

void FeClusterReactionNetwork::reinitializeNetwork() {

	// Reset the Ids
	// std::for_each is guaranteed to visit reactants in order for C++11.
	int id = 0;
	std::for_each(allReactants.begin(), allReactants.end(),
			[&id](IReactant& currReactant) {
				id++;
				currReactant.setId(id);
				currReactant.setMomentId(id, 0);
				currReactant.setMomentId(id, 1);
			});

	// Get all the super clusters and loop on them
	// Have to use allReactants again to be sure the ordering is the same across plateforms
	std::for_each(allReactants.begin(), allReactants.end(),
			[&id, this](IReactant& currReactant) {
				if (currReactant.getType() == ReactantType::FeSuper) {

					auto& currCluster = static_cast<FeSuperCluster&>(currReactant);

					id++;
					currCluster.setMomentId(id, 0);
					id++;
					currCluster.setMomentId(id, 1);

					// Update the HeV size
					auto const& heBounds = currCluster.getHeBounds();
					auto const& vBounds = currCluster.getVBounds();
					IReactant::SizeType clusterSize = (*(heBounds.end()) - 1)
					+ (*(vBounds.end()) - 1);
					if (clusterSize > maxClusterSizeMap[ReactantType::HeV]) {
						maxClusterSizeMap[ReactantType::HeV] = clusterSize;
					}
				}
			});

			return;
		}

		void FeClusterReactionNetwork::reinitializeConnectivities() {

			// Reset connectivities of each reactant.
			std::for_each(allReactants.begin(), allReactants.end(),
					[](IReactant& currReactant) {
						currReactant.resetConnectivities();
					});

			return;
		}

		void FeClusterReactionNetwork::updateConcentrationsFromArray(
				double * concentrations) {

			// Set the concentration on each reactant.
			std::for_each(allReactants.begin(), allReactants.end(),
					[&concentrations](IReactant& currReactant) {
						auto id = currReactant.getId() - 1;
						currReactant.setConcentration(concentrations[id]);
					});

			// Set the moments
			auto const& superTypeMap = getAll(ReactantType::FeSuper);
			std::for_each(superTypeMap.begin(), superTypeMap.end(),
					[&concentrations](const ReactantMap::value_type& currMapItem) {

						auto& cluster = static_cast<FeSuperCluster&>(*(currMapItem.second));

						cluster.setZerothMoment(concentrations[cluster.getId() - 1]);
						cluster.setHeMoment(concentrations[cluster.getMomentId(0) - 1]);
						cluster.setVMoment(concentrations[cluster.getMomentId(1) - 1]);
					});

			return;
		}

		std::vector<std::vector<int> > FeClusterReactionNetwork::getCompositionList() const {
			// Create the list that will be returned
			std::vector<std::vector<int> > compList;

			// Loop on all the reactants
			std::for_each(allReactants.begin(), allReactants.end(),
					[&compList](IReactant& currReactant) {
						// Get the composition
						auto comp = currReactant.getComposition();
						std::vector <int> compVec;
						compVec.push_back(comp[toCompIdx(Species::He)]);
						compVec.push_back(comp[toCompIdx(Species::V)]);
						compVec.push_back(comp[toCompIdx(Species::I)]);

						// Save the composition in the list
						compList.push_back(compVec);
					});

			return compList;
		}

		void FeClusterReactionNetwork::getDiagonalFill(SparseFillMap& fillMap) {
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
							if (connectivity[j] == 1) {
								fillMap[id].emplace_back(j);
								columnIds.push_back(j);
							}
						}
						// Update the map
						dFillMap[id] = columnIds;
					});

			// Get the connectivity for each moment
			for (auto const& superMapItem : getAll(ReactantType::FeSuper)) {

				auto const& reactant =
						static_cast<FeSuperCluster&>(*(superMapItem.second));

				// Get the reactant and its connectivity
				auto const& connectivity = reactant.getConnectivity();
				auto connectivityLength = connectivity.size();
				// Get the helium moment id so that the connectivity can be lined up in
				// the proper column
				auto id = reactant.getMomentId(0) - 1;

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

				// Get the vacancy moment id so that the connectivity can be lined up in
				// the proper column
				id = reactant.getMomentId(1) - 1;

				// Add it to the diagonal fill block
				for (int j = 0; j < connectivityLength; j++) {
					if (connectivity[j] == 1) {
						fillMap[id].emplace_back(j);
					}
				}
				// Update the map
				dFillMap[id] = columnIds;
			}

			return;
		}

		double FeClusterReactionNetwork::getTotalAtomConcentration(int i) {
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
				heliumConc += cluster.getConcentration()
						* comp[toCompIdx(Species::He)];
			}

			// Sum over all super clusters.
			for (auto const& currMapItem : getAll(ReactantType::FeSuper)) {

				// Get the cluster
				auto const& cluster =
						static_cast<FeSuperCluster&>(*(currMapItem.second));

				// Add its total helium concentration helium concentration
				heliumConc += cluster.getTotalHeliumConcentration();
			}

			return heliumConc;
		}

		double FeClusterReactionNetwork::getTotalTrappedAtomConcentration(
				int i) {
			// Initial declarations
			double heliumConc = 0.0;

			// Sum over all HeV clusters.
			for (auto const& currMapItem : getAll(ReactantType::HeV)) {
				// Get the cluster and its composition
				auto const& cluster = *(currMapItem.second);
				auto& comp = cluster.getComposition();

				// Add the concentration times the He content to the total helium concentration
				heliumConc += cluster.getConcentration()
						* comp[toCompIdx(Species::He)];
			}

			// Sum over all super clusters.
			for (auto const& currMapItem : getAll(ReactantType::FeSuper)) {
				// Get the cluster
				auto const& cluster =
						static_cast<FeSuperCluster&>(*(currMapItem.second));

				// Add its total helium concentration
				heliumConc += cluster.getTotalHeliumConcentration();
			}

			return heliumConc;
		}

		double FeClusterReactionNetwork::getTotalVConcentration() {
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
				vConc += cluster.getConcentration()
						* comp[toCompIdx(Species::V)];
			}

			// Sum over all super clusters
			for (auto const& currMapItem : getAll(ReactantType::FeSuper)) {
				// Get the cluster
				auto const& cluster =
						static_cast<FeSuperCluster&>(*(currMapItem.second));

				// Add its total vacancy concentration
				vConc += cluster.getTotalVacancyConcentration();
			}

			return vConc;
		}

		double FeClusterReactionNetwork::getTotalIConcentration() {
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

		void FeClusterReactionNetwork::computeAllFluxes(
				double *updatedConcOffset, int i) {

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
			for (auto const& currMapItem : getAll(ReactantType::FeSuper)) {

				auto const& superCluster =
						static_cast<FeSuperCluster&>(*(currMapItem.second));

				// Compute the helium moment flux
				auto flux = superCluster.getHeMomentFlux();
				// Update the concentration of the cluster
				auto reactantIndex = superCluster.getMomentId(0) - 1;
				updatedConcOffset[reactantIndex] += flux;

				// Compute the vacancy moment flux
				flux = superCluster.getVMomentFlux();
				// Update the concentration of the cluster
				reactantIndex = superCluster.getMomentId(1) - 1;
				updatedConcOffset[reactantIndex] += flux;
			}

			return;
		}

		void FeClusterReactionNetwork::computeAllPartials(
				const std::vector<size_t>& startingIdx,
				const std::vector<int>& indices, std::vector<double>& vals,
				int i) const {
			// Initial declarations
			const int dof = getDOF();
			std::vector<double> clusterPartials(dof, 0.0);

			// Get the super clusters
			auto const& superClusters = getAll(ReactantType::FeSuper);

			// Make a vector of types for the non super clusters
			std::vector<ReactantType> typeVec { ReactantType::He,
					ReactantType::V, ReactantType::I, ReactantType::HeV };
			// Loop on it
			for (auto tvIter = typeVec.begin(); tvIter != typeVec.end();
					++tvIter) {

				auto currType = *tvIter;

				// Consider all reactants of the current type.
				auto const& currTypeReactantMap = getAll(currType);

				// Update the column in the Jacobian that represents each normal reactant
				for (auto const& currMapItem : currTypeReactantMap) {

					auto const& reactant =
							static_cast<FeCluster&>(*(currMapItem.second));

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
						vals[myStartingIdx + j] =
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
						static_cast<FeSuperCluster&>(*(currMapItem.second));

				{
					// Get the super cluster index
					auto reactantIndex = reactant.getId() - 1;

					// Get the partial derivatives
					reactant.getPartialDerivatives(clusterPartials, i);

					// Get the list of column ids from the map
					auto const& pdColIdsVector = dFillMap.at(reactantIndex);

					// Loop over the list of column ids
					auto myStartingIdx = startingIdx[reactantIndex];
					for (int j = 0; j < pdColIdsVector.size(); j++) {
						// Get the partial derivative from the array of all of the partials
						vals[myStartingIdx + j] =
								clusterPartials[pdColIdsVector[j]];

						// Reset the cluster partial value to zero. This is much faster
						// than using memset.
						clusterPartials[pdColIdsVector[j]] = 0.0;
					}
				}

				{
					// Get the helium moment index
					auto reactantIndex = reactant.getMomentId(0) - 1;

					// Get the partial derivatives
					reactant.getHeMomentPartialDerivatives(clusterPartials);
					// Get the list of column ids from the map
					auto const& pdColIdsVector = dFillMap.at(reactantIndex);

					// Loop over the list of column ids
					auto myStartingIdx = startingIdx[reactantIndex];
					for (int j = 0; j < pdColIdsVector.size(); j++) {
						// Get the partial derivative from the array of all of the partials
						vals[myStartingIdx + j] =
								clusterPartials[pdColIdsVector[j]];

						// Reset the cluster partial value to zero. This is much faster
						// than using memset.
						clusterPartials[pdColIdsVector[j]] = 0.0;
					}
				}

				{
					// Get the vacancy moment index
					auto reactantIndex = reactant.getMomentId(1) - 1;

					// Get the partial derivatives
					reactant.getVMomentPartialDerivatives(clusterPartials);
					// Get the list of column ids from the map
					auto const& pdColIdsVector = dFillMap.at(reactantIndex);

					// Loop over the list of column ids
					auto myStartingIdx = startingIdx[reactantIndex];
					for (int j = 0; j < pdColIdsVector.size(); j++) {
						// Get the partial derivative from the array of all of the partials
						vals[myStartingIdx + j] =
								clusterPartials[pdColIdsVector[j]];

						// Reset the cluster partial value to zero. This is much faster
						// than using memset.
						clusterPartials[pdColIdsVector[j]] = 0.0;
					}
				}
			}

			return;
		}

		double FeClusterReactionNetwork::computeBindingEnergy(
				const DissociationReaction& reaction) const {

			double bindingEnergy = 5.0;
			if (reaction.dissociating.getType() == ReactantType::He
					&& reaction.first.getType() == ReactantType::He) {
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
			if ((reaction.dissociating.getType() == ReactantType::HeV)
					&& (reaction.first.getType() == ReactantType::V
							|| reaction.second.getType() == ReactantType::V)) {
				auto& comp = reaction.dissociating.getComposition();
				bindingEnergy =
						1.73
								- 2.59
										* (pow(
												(double) comp[toCompIdx(
														Species::V)], 2.0 / 3.0)
												- pow(
														(double) comp[toCompIdx(
																Species::V)]
																- 1.0,
														2.0 / 3.0))
								+ 2.5
										* log(
												1.0
														+ ((double) comp[toCompIdx(
																Species::He)]
																/ (double) comp[toCompIdx(
																		Species::V)]));
			}
			if (reaction.dissociating.getType() == ReactantType::FeSuper
					&& (reaction.first.getType() == ReactantType::V
							|| reaction.second.getType() == ReactantType::V)) {
				auto& comp = reaction.dissociating.getComposition();
				double numV = (double) comp[toCompIdx(Species::V)];
				double numHe = (double) comp[toCompIdx(Species::He)];
				bindingEnergy = 1.73
						- 2.59
								* (pow(numV, 2.0 / 3.0)
										- pow(numV - 1.0, 2.0 / 3.0))
						+ 2.5 * log(1.0 + (numHe / numV));
			}
			if (reaction.first.getType() == ReactantType::I
					|| reaction.second.getType() == ReactantType::I) {
				if (reaction.dissociating.getType() == ReactantType::HeV) {
					auto& comp = reaction.dissociating.getComposition();
					bindingEnergy =
							4.88
									+ 2.59
											* (pow(
													(double) comp[toCompIdx(
															Species::V)],
													2.0 / 3.0)
													- pow(
															(double) comp[toCompIdx(
																	Species::V)]
																	- 1.0,
															2.0 / 3.0))
									- 2.5
											* log(
													1.0
															+ ((double) comp[toCompIdx(
																	Species::He)]
																	/ (double) comp[toCompIdx(
																			Species::V)]));
				} else if (reaction.dissociating.getType()
						== ReactantType::FeSuper) {
					auto& comp = reaction.dissociating.getComposition();
					double numV = (double) comp[toCompIdx(Species::V)];
					double numHe = (double) comp[toCompIdx(Species::He)];
					bindingEnergy = 4.88
							+ 2.59
									* (pow(numV, 2.0 / 3.0)
											- pow(numV - 1.0, 2.0 / 3.0))
							- 2.5 * log(1.0 + (numHe / numV));
				} else if (reaction.dissociating.getType()
						== ReactantType::He) {
					int size = reaction.dissociating.getSize();
					switch (size) {
					case 1:
						bindingEnergy = 4.31;
						break;
					case 2:
						bindingEnergy = 2.90;
						break;
					case 3:
						bindingEnergy = 2.02;
						break;
					case 4:
						bindingEnergy = 1.09;
						break;
					case 5:
						bindingEnergy = 0.58;
						break;
					case 6:
						bindingEnergy = 0.13;
						break;
					case 7:
						bindingEnergy = -0.25;
						break;
					case 8:
						bindingEnergy = -0.59;
						break;
					default:
						break;
					}
				}

			}

//	if (bindingEnergy < -5.0)
//	std::cout << "dissociation: " << reaction.dissociating.getName() << " -> "
//			<< reaction.first.getName() << " + "
//			<< reaction.second.getName() << " : " << bindingEnergy
//			<< std::endl;

			return max(bindingEnergy, -5.0);
		}

		IReactant * FeClusterReactionNetwork::getSuperFromComp(
				IReactant::SizeType nHe, IReactant::SizeType nV) {

			// Requests for finding a particular supercluster have high locality.
			// See if the last supercluster we were asked to find is the right
			// one for this request.
			static IReactant* lastRet = nullptr;
			if (lastRet
					and static_cast<FeSuperCluster*>(lastRet)->isIn(nHe, nV)) {
				return lastRet;
			}

			// We didn't find the last supercluster in our cache, so do a full lookup.
			IReactant* ret = nullptr;

			auto heBaseIdx = findBoundsIntervalBaseIdx(nHe);
			auto vBaseIdx = findBoundsIntervalBaseIdx(nV);

			if ((heBaseIdx != std::numeric_limits<std::size_t>::max())
					and (vBaseIdx != std::numeric_limits<std::size_t>::max())) {

				auto& superIter = superClusterLookupMap[heBaseIdx][vBaseIdx];
				if (superIter
						!= clusterTypeMap.at(ReactantType::FeSuper).end()) {
					ret = superIter->second.get();
					assert(static_cast<FeSuperCluster*>(ret)->isIn(nHe, nV));
					lastRet = ret;
				}
			}

			return ret;
		}

		} // namespace xolotlCore

