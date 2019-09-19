#include "AlloyClusterReactionNetwork.h"
#include "AlloySuperCluster.h"
#include <xolotlPerf.h>
#include <sstream>
#include <Constants.h>
#include "AlloyCases.h"

using namespace xolotlCore;

AlloyClusterReactionNetwork::AlloyClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork( { ReactantType::V, ReactantType::I,
				ReactantType::Frank, ReactantType::Void, ReactantType::Perfect,
				ReactantType::Faulted, ReactantType::VoidSuper,
				ReactantType::FaultedSuper, ReactantType::PerfectSuper,
				ReactantType::FrankSuper }, registry) {

	// Initialize default properties
	dissociationsEnabled = true;

	return;
}

double AlloyClusterReactionNetwork::calculateReactionRateConstant(
		const ProductionReaction& reaction, int i) const {
	// Get the reaction radii
	double r_first = reaction.first.getReactionRadius();
	double r_second = reaction.second.getReactionRadius();

	// Get the diffusion coefficients
	double firstDiffusion = reaction.first.getDiffusionCoefficient(i);
	double secondDiffusion = reaction.second.getDiffusionCoefficient(i);
	
	// Initialize rate constant and parameters to calculate k+ for loop/sphere case
	double k_plus, zl, intBias, p = 0.0;
	double zs = 4.0 * xolotlCore::pi * 
		(r_first + r_second + xolotlCore::alloyCoreRadius);

	// Calculate and return
	// Formula for reaction rate constant depends on type of reactants
	
	// k+ for reaction between two spheres
	if (reaction.first.isSphere() && reaction.second.isSphere()) {
		k_plus = 4.0 * xolotlCore::pi
			* (r_first + r_second + xolotlCore::alloyCoreRadius)
			* (firstDiffusion + secondDiffusion);
	}

        // k+ for reaction between loop and sphere or between two loops
        else {
                if (r_first >= r_second) {
                        p = 1.0 / (1.0 + pow(r_first /
                        (r_second + xolotlCore::alloyCoreRadius), 2.0));
                        zl = 4.0 * pow(xolotlCore::pi, 2.0) * r_first
                               / log(1.0 + 8.0 * r_first
                               / (r_second + xolotlCore::alloyCoreRadius));
                        }
                else {
                        p = 1.0 / (1.0 + pow(r_second /
                               (r_first + xolotlCore::alloyCoreRadius), 2.0));
                        zl = 4.0 * pow(xolotlCore::pi, 2.0) * r_second
                                / log(1.0 + 8.0 * r_second
                                / (r_first + xolotlCore::alloyCoreRadius));
                        }
                k_plus = (firstDiffusion + secondDiffusion) * (p * zs
                                + (1.0 - p) / zl);
        }
	
	// Interstitial bias for interactions with loops
	if ((reaction.first.getType() == ReactantType::I and not reaction.second.isSphere())
           || (reaction.second.getType() == ReactantType::I and not reaction.first.isSphere())) {
		intBias = 1.2;
	}
	
	else {
		intBias = 1.0;
	}
	
	k_plus *= intBias;
	
	return k_plus;
}

double AlloyClusterReactionNetwork::calculateDissociationConstant(
		const DissociationReaction& reaction, int i) {

	// If the dissociations are not allowed
	if (!dissociationsEnabled)
		return 0.0;

	// Compute the atomic volume (included the fact that there are 4 atoms per cell)
	double atomicVolume = 0.25 * xolotlCore::alloyLatticeConstant
			* xolotlCore::alloyLatticeConstant
			* xolotlCore::alloyLatticeConstant;

	// Get the rate constant from the reverse reaction
	double kPlus = reaction.reverseReaction->kConstant[i];

	// Calculate Binding Energy
	double bindingEnergy = computeBindingEnergy(reaction);

	// Correct smallest faulted loop binding energy
	int minFaultedSize = 7;
	if (reaction.dissociating.getType() == ReactantType::Faulted
			&& reaction.dissociating.getSize() == minFaultedSize) {
		bindingEnergy = 1.9
				- 3.2
						* (pow(double(minFaultedSize), 2.0 / 3.0)
								- pow(double(minFaultedSize - 1), 2.0 / 3.0));
	}

	double k_minus_exp = exp(
			-1.0 * bindingEnergy / (xolotlCore::kBoltzmann * temperature));	// We can use the network temperature
	// because this method is called only
	// when the temperature is updated
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;

//	std::cout << reaction.first.getName() << " + " << reaction.second.getName()
//			<< " <- " << reaction.dissociating.getName() << " : "
//			<< k_minus << " " << bindingEnergy << std::endl;

	return k_minus;
}

double AlloyClusterReactionNetwork::computeBindingEnergy(
		const DissociationReaction& reaction) const {
//	double firstEnergy = reaction.first.getFormationEnergy(), secondEnergy =
//			reaction.second.getFormationEnergy(), dissoEnergy =
//			reaction.dissociating.getFormationEnergy();
//
//	auto& dissoCluster = static_cast<AlloyCluster&>(reaction.dissociating);
//	if (dissoCluster.isSuper()) {
//		dissoEnergy = getFormationEnergy(dissoCluster.getType(),
//				dissoCluster.getSize());
//		if (reaction.first.getSize() == 1)
//			secondEnergy = getFormationEnergy(dissoCluster.getType(),
//					dissoCluster.getSize() - 1);
//		else
//			firstEnergy = getFormationEnergy(dissoCluster.getType(),
//					dissoCluster.getSize() - 1);
//	}
//
	double bindingEnergy = 0.0;

	if (reaction.dissociating.getType() == ReactantType::Void
		|| reaction.dissociating.getType() == ReactantType::VoidSuper) {
			if (reaction.first.getType() == ReactantType::I
					|| reaction.second.getType() == ReactantType::I) {
				double n = reaction.dissociating.getSize();
				bindingEnergy = 3.5
						- 3.45 * (pow(n + 1.0, 2.0 / 3.0) - pow(n, 2.0 / 3.0));
			} else if (reaction.first.getType() == ReactantType::V
					|| reaction.second.getType() == ReactantType::V) {
				double n = reaction.dissociating.getSize();
				bindingEnergy = 1.9
						- 3.1 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
			}
		} else if (reaction.dissociating.getType() == ReactantType::Faulted
				|| reaction.dissociating.getType()
						== ReactantType::FaultedSuper) {
			if (reaction.first.getType() == ReactantType::V
					|| reaction.second.getType() == ReactantType::V) {
				double n = reaction.dissociating.getSize();
				bindingEnergy = 1.9
						- 3.2 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
			}
		} else if (reaction.dissociating.getType() == ReactantType::V) {
			if (reaction.first.getType() == ReactantType::V 
					|| reaction.second.getType() == ReactantType::V) {
				double n = reaction.dissociating.getSize();
				bindingEnergy = 1.9 
						- 3.1 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
			}
		} else if (reaction.dissociating.getType() == ReactantType::I) {
			if (reaction.first.getType() == ReactantType::I
					|| reaction.second.getType() == ReactantType::I) {
				double n = reaction.dissociating.getSize();
				bindingEnergy = 3.5
						- 2.5 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
			}
		}

//	std::cout << reaction.dissociating.getName() << " -> "
//			<< reaction.first.getName() << " + " << reaction.second.getName()
//			<< " : " << bindingEnergy << " " << dissoEnergy << " "
//			<< firstEnergy << " " << secondEnergy << std::endl;
//
	return bindingEnergy;
}

int AlloyClusterReactionNetwork::typeSwitch(ReactantType const typeName) const {
	if (typeName == ReactantType::V || typeName == ReactantType::Void
			|| typeName == ReactantType::Faulted
			|| typeName == ReactantType::VoidSuper
			|| typeName == ReactantType::FaultedSuper)
		return -1;
	else
		return 1;
}

double AlloyClusterReactionNetwork::getReactionRadius(
		ReactantType const typeName, int size) const {
	if (typeName == ReactantType::Perfect || typeName == ReactantType::PerfectSuper)
		return sqrt((double(size) * 0.25 * xolotlCore::alloyLatticeConstant 
			* xolotlCore::alloyLatticeConstant * xolotlCore::alloyLatticeConstant) 
			/ (xolotlCore::pi * xolotlCore::alloyLatticeConstant
			  * xolotlCore::perfectBurgers));
	if (typeName == ReactantType::Faulted || typeName == ReactantType::FaultedSuper)
		return sqrt((double(size) * 0.25 * xolotlCore::alloyLatticeConstant 
			* xolotlCore::alloyLatticeConstant * xolotlCore::alloyLatticeConstant) 
			/ (xolotlCore::pi * xolotlCore::alloyLatticeConstant
			  * xolotlCore::faultedBurgers));
	if (typeName == ReactantType::Frank || typeName == ReactantType::FrankSuper)
		return sqrt((double(size) * 0.25 * xolotlCore::alloyLatticeConstant 
			* xolotlCore::alloyLatticeConstant * xolotlCore::alloyLatticeConstant) 
			/ (xolotlCore::pi * xolotlCore::alloyLatticeConstant
			  * xolotlCore::frankBurgers));
	if (typeName == ReactantType::V || typeName == ReactantType::Void
			|| typeName == ReactantType::I
			|| typeName == ReactantType::VoidSuper)
		return pow((0.75 * 0.25 * xolotlCore::alloyLatticeConstant * 
			xolotlCore::alloyLatticeConstant * xolotlCore::alloyLatticeConstant 
			* double(size)) / xolotlCore::pi, 1.0 / 3.0);

	return 0.0;
}

double AlloyClusterReactionNetwork::getFormationEnergy(
		ReactantType const typeName, int size) const {
	if (typeName == ReactantType::Perfect || typeName == ReactantType::Frank
			|| typeName == ReactantType::PerfectSuper
			|| typeName == ReactantType::FrankSuper)
		return 3.4 + 2.0 * (pow(double(size), 2.0 / 3.0) - 1.0);
	if (typeName == ReactantType::Faulted
			|| typeName == ReactantType::FaultedSuper)
		return 1.9 + 2.0 * (pow(double(size), 2.0 / 3.0) - 1.0);
	if (typeName == ReactantType::Void || typeName == ReactantType::VoidSuper)
		return 1.9 + 3.4 * (pow(double(size), 2.0 / 3.0) - 1.0);
	if (typeName == ReactantType::V)
		return 1.9 + 3.4 * (pow(double(size), 2.0 / 3.0) - 1.0);
	if (typeName == ReactantType::I)
		return 3.4 + 3.5 * (pow(double(size), 2.0 / 3.0) - 1.0);
	return 0.0;
}

void AlloyClusterReactionNetwork::createReactionConnectivity() {
	// Initial declarations
	int firstSize = 0, secondSize = 0, productSize = 0;

	// Get forward reaction network
	auto forwardReactions = getForwardReactions("default");

	// Loop over all reactions
	for (auto & forwardReaction : forwardReactions) {
		// Get all reactants for given reaction
		auto& allReactants1 = getAll(forwardReaction.getFirstReactant());
		auto& allReactants2 = getAll(forwardReaction.getSecondReactant());
		// Loop over all individual reactants
		for (auto & reactant1 : allReactants1) {
			auto& cluster1 = static_cast<AlloyCluster&>(*(reactant1.second));
			for (auto & reactant2 : allReactants2) {
				auto& cluster2 = static_cast<AlloyCluster&>(*(reactant2.second));
				// Skip repeating reactions
				if ((cluster1.getType() == cluster2.getType())
						&& (cluster2.getSize() > cluster1.getSize()))
					continue;
				// Skip if both are immobile
				if ((cluster1.getDiffusionFactor() == 0.0)
						&& (cluster2.getDiffusionFactor() == 0.0))
					continue;

				// If both clusters are not super
				if (!cluster1.isSuper() && !cluster2.isSuper()) {
					// Look for recombination
					int size1 = cluster1.getSize()
							* typeSwitch(cluster1.getType());
					int size2 = cluster2.getSize()
							* typeSwitch(cluster2.getType());
					int productSize = size1 + size2;
					if (productSize == 0) {
						// Create the reaction
						std::unique_ptr<ProductionReaction> reaction(
								new ProductionReaction(cluster1, cluster2));
						auto& prref = add(std::move(reaction));
						// Tell the reactants that they are in this reaction
						(cluster1).participateIn(prref);
						(cluster2).participateIn(prref);

//								std::cout << cluster1.getName() << " + "
//										<< cluster2.getName() << " -> "
//										<< "recombined" << std::endl;
					}
				}

				// Get list of accepted products
				auto products = forwardReaction.getProducts();
				// Loop over all accepted products
				for (auto & productName : products) {
					// Get all the clusters of this type
					auto& allProducts = getAll(productName);
					// Loop over all individual products
					for (auto & product : allProducts) {
						auto& prodCluster =
								static_cast<AlloyCluster&>(*(product.second));
						// Check if the reaction is possible
						if (checkOverlap(cluster1, cluster2, prodCluster)) {
							// Create the reaction
							std::unique_ptr<ProductionReaction> reaction(
									new ProductionReaction(cluster1, cluster2));
							auto& prref = add(std::move(reaction));
							// Tell the reactants that they are in this reaction
							(cluster1).participateIn(prref, prodCluster);
							(cluster2).participateIn(prref, prodCluster);
							prodCluster.resultFrom(prref, prodCluster);
						}
					}
				}
			}
		}
	}

	auto backwardReactions = getBackwardReactions("default");
	for (auto & backwardReaction : backwardReactions) {
		auto monomerName = backwardReaction.getMonomer();
		auto& monomer = static_cast<AlloyCluster&>(*get(toSpecies(monomerName),
				1));
		auto parentName = backwardReaction.getParent();
		auto& parents = getAll(parentName);
		for (auto & parent : parents) {
			auto& parentCluster = static_cast<AlloyCluster&>(*(parent.second));
			auto productNames = backwardReaction.getProducts();
			// Loop over all accepted products
			for (auto & productName : productNames) {
				// Get all the clusters of this type
				auto& allProducts = getAll(productName);
				// Loop over all individual products
				for (auto & product : allProducts) {
					auto& prodCluster =
							static_cast<AlloyCluster&>(*(product.second));
					if (checkOverlap(monomer, prodCluster, parentCluster)) {
						// Create the dissociation
						std::unique_ptr<DissociationReaction> dissoReaction(
								new DissociationReaction(parentCluster, monomer,
										prodCluster));
						// Create the corresponding reaction
						std::unique_ptr<ProductionReaction> reaction(
								new ProductionReaction(monomer, prodCluster));
						auto& prref = add(std::move(reaction));
						// Set it in the dissociation
						dissoReaction->reverseReaction = &prref;
						auto& drref = add(std::move(dissoReaction));
						// Tell the reactants that they are in this reaction
						(monomer).participateIn(drref, parentCluster);
						(prodCluster).participateIn(drref, parentCluster);
						parentCluster.emitFrom(drref, parentCluster);
					}
				}
			}
		}
	}

	return;
}

void AlloyClusterReactionNetwork::checkDissociationConnectivity(
		IReactant * emittingReactant,
		std::shared_ptr<ProductionReaction> reaction) {

	return;
}

void AlloyClusterReactionNetwork::setTemperature(double temp, int i) {
	ReactionNetwork::setTemperature(temp, i);

	computeRateConstants(i);

	return;
}

void AlloyClusterReactionNetwork::reinitializeNetwork() {
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

				if (currReactant.getType() == ReactantType::VoidSuper
						|| currReactant.getType() == ReactantType::FaultedSuper
						|| currReactant.getType() == ReactantType::FrankSuper
						|| currReactant.getType() == ReactantType::PerfectSuper) {
					auto& currCluster = static_cast<AlloySuperCluster&>(currReactant);
					id++;
					currCluster.setMomentId(id);

					// Update the size
					IReactant::SizeType clusterSize = (double)currCluster.getSize()
					+ (double)(currCluster.getSectionWidth() - 1) / 2.0;
					if (currReactant.getType() == ReactantType::VoidSuper
							&& clusterSize > maxClusterSizeMap[ReactantType::Void]) {
						maxClusterSizeMap[ReactantType::Void] = clusterSize;
					}
					if (currReactant.getType() == ReactantType::FaultedSuper
							&& clusterSize > maxClusterSizeMap[ReactantType::Faulted]) {
						maxClusterSizeMap[ReactantType::Faulted] = clusterSize;
					}
					if (currReactant.getType() == ReactantType::PerfectSuper
							&& clusterSize > maxClusterSizeMap[ReactantType::Perfect]) {
						maxClusterSizeMap[ReactantType::Perfect] = clusterSize;
					}
					if (currReactant.getType() == ReactantType::FrankSuper
							&& clusterSize > maxClusterSizeMap[ReactantType::Frank]) {
						maxClusterSizeMap[ReactantType::Frank] = clusterSize;
					}
				}
			});

	return;
}

void AlloyClusterReactionNetwork::reinitializeConnectivities() {
	// Loop on all the reactants to reset their connectivities
	std::for_each(allReactants.begin(), allReactants.end(),
			[](IReactant& currReactant) {
				currReactant.resetConnectivities();
			});

	return;
}

void AlloyClusterReactionNetwork::updateConcentrationsFromArray(
		double * concentrations) {

	// Set the concentration on each reactant.
	std::for_each(allReactants.begin(), allReactants.end(),
			[&concentrations](IReactant& currReactant) {
				auto id = currReactant.getId() - 1;
				currReactant.setConcentration(concentrations[id]);
				// Set the moments
				if (currReactant.getType() == ReactantType::VoidSuper
						|| currReactant.getType() == ReactantType::FaultedSuper
						|| currReactant.getType() == ReactantType::FrankSuper
						|| currReactant.getType() == ReactantType::PerfectSuper) {
					auto& currCluster = static_cast<AlloySuperCluster&>(currReactant);
					currCluster.setZerothMoment(concentrations[id]);
					currCluster.setMoment(concentrations[currCluster.getMomentId() - 1]);
				}
			});

	return;
}

std::vector<std::vector<int> > AlloyClusterReactionNetwork::getCompositionList() const {
	// Create the list that will be returned
	std::vector<std::vector<int> > compList;

	// Loop on all the reactants
	std::for_each(allReactants.begin(), allReactants.end(),
			[&compList](IReactant& currReactant) {
				// Get the composition
				auto comp = currReactant.getComposition();
				std::vector <int> compVec;
				compVec.push_back(comp[toCompIdx(Species::V)]);
				compVec.push_back(comp[toCompIdx(Species::I)]);
				compVec.push_back(comp[toCompIdx(Species::Void)]);
				compVec.push_back(comp[toCompIdx(Species::Perfect)]);
				compVec.push_back(comp[toCompIdx(Species::Faulted)]);
				compVec.push_back(comp[toCompIdx(Species::Frank)]);

				// Save the composition in the list
				compList.push_back(compVec);
			});

	return compList;
}

void AlloyClusterReactionNetwork::getDiagonalFill(SparseFillMap& fillMap) {

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
	// Make a vector of types for the super clusters
	std::vector<ReactantType> typeVec { ReactantType::PerfectSuper,
			ReactantType::FaultedSuper, ReactantType::FrankSuper,
			ReactantType::VoidSuper };
	// Loop on it
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {

		auto currType = *tvIter;

		// Consider all reactants of the current type.
		auto const& currTypeReactantMap = getAll(currType);

		// Update the column in the Jacobian that represents each normal reactant
		for (auto const& currMapItem : currTypeReactantMap) {

			auto const& reactant =
					static_cast<AlloySuperCluster&>(*(currMapItem.second));

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
	}

	return;
}

void AlloyClusterReactionNetwork::computeAllFluxes(double *updatedConcOffset,
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
	// Make a vector of types for the super clusters
	std::vector<ReactantType> typeVec { ReactantType::PerfectSuper,
			ReactantType::FaultedSuper, ReactantType::FrankSuper,
			ReactantType::VoidSuper };
	// Loop on it
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {

		auto currType = *tvIter;

		// Consider all reactants of the current type.
		auto const& currTypeReactantMap = getAll(currType);

		// Update the column in the Jacobian that represents each normal reactant
		for (auto const& currMapItem : currTypeReactantMap) {

			auto& superCluster =
					static_cast<AlloySuperCluster&>(*(currMapItem.second));

			// Compute the xenon moment flux
			auto flux = superCluster.getMomentFlux();
			// Update the concentration of the cluster
			auto reactantIndex = superCluster.getMomentId() - 1;
			updatedConcOffset[reactantIndex] += flux;
		}
	}

	return;
}

void AlloyClusterReactionNetwork::computeAllPartials(
		const std::vector<size_t>& startingIdx, const std::vector<int>& indices,
		std::vector<double>& vals, int i) const {

	// Initial declarations
	const int dof = getDOF();
	std::vector<double> clusterPartials(dof, 0.0);

	// Make a vector of types for the non super clusters
	std::vector<ReactantType> typeVec { ReactantType::V, ReactantType::I,
			ReactantType::Void, ReactantType::Perfect, ReactantType::Faulted,
			ReactantType::Frank };
	// Loop on it
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {

		auto currType = *tvIter;

		// Consider all reactants of the current type.
		auto const& currTypeReactantMap = getAll(currType);

		// Update the column in the Jacobian that represents each normal reactant
		for (auto const& currMapItem : currTypeReactantMap) {

			auto const& reactant =
					static_cast<AlloyCluster&>(*(currMapItem.second));

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
	}

	// Update the column in the Jacobian that represents the moment for the super clusters
	// Make a vector of types for the super clusters
	std::vector<ReactantType> superTypeVec { ReactantType::PerfectSuper,
			ReactantType::FaultedSuper, ReactantType::FrankSuper,
			ReactantType::VoidSuper };
	// Loop on it
	for (auto tvIter = superTypeVec.begin(); tvIter != superTypeVec.end();
			++tvIter) {

		auto currType = *tvIter;

		// Consider all reactants of the current type.
		auto const& currTypeReactantMap = getAll(currType);

		// Update the column in the Jacobian that represents each normal reactant
		for (auto const& currMapItem : currTypeReactantMap) {

			auto const& reactant =
					static_cast<AlloySuperCluster&>(*(currMapItem.second));

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
					vals[myStartingIdx + j] =
							clusterPartials[pdColIdsVector[j]];

					// Reset the cluster partial value to zero. This is much faster
					// than using memset.
					clusterPartials[pdColIdsVector[j]] = 0.0;
				}
			}
			{
				// Get the moment index
				auto reactantIndex = reactant.getMomentId() - 1;

				// Get the partial derivatives
				reactant.getMomentPartialDerivatives(clusterPartials);
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
	}

	return;
}
