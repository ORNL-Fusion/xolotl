// Includes
#include "NESuperCluster.h"
#include "NEClusterReactionNetwork.h"
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

/**
 * The xenon moment partials.
 */
std::vector<double> momentPartials;

NESuperCluster::NESuperCluster(double num, int nTot, int width, double radius,
		double energy, IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		NECluster(_network, registry, buildName(num)), numXe(num), nTot(nTot), l0(
				0.0), l1(0.0), dispersion(0.0), momentFlux(0.0) {
	// Set the cluster size
	size = (int) numXe;

	// Update the composition map
	composition[toCompIdx(Species::Xe)] = (int) (numXe * (double) nTot);

	// Set the width
	sectionWidth = width;

	// Set the reaction radius and formation energy
	reactionRadius = radius;
	formationEnergy = energy;

	// Set the diffusion factor and the migration energy
	migrationEnergy = std::numeric_limits<double>::infinity();
	diffusionFactor = 0.0;

	// Set the typename appropriately
	type = ReactantType::NESuper;

	return;
}

void NESuperCluster::resultFrom(ProductionReaction& reaction, double *coef) {

	// Create a new SuperClusterProductionPair
	SuperClusterProductionPair superPair(
			&static_cast<NECluster&>(reaction.first),
			&static_cast<NECluster&>(reaction.second), &reaction);
	// Update the coeficients
	superPair.a000 = coef[0];
	superPair.a001 = coef[1];
	superPair.a100 = coef[2];
	superPair.a101 = coef[3];
	superPair.a010 = coef[4];
	superPair.a011 = coef[5];
	superPair.a110 = coef[6];
	superPair.a111 = coef[7];

	// Add it to the list
	effReactingList.push_front(superPair);

	return;
}

void NESuperCluster::participateIn(ProductionReaction& reaction, double *coef) {
	// Look for the other cluster
	auto& otherCluster = static_cast<NECluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Create a new SuperClusterProductionPair with NULL as the second cluster because
	// we do not need it
	SuperClusterProductionPair superPair(&otherCluster, nullptr, &reaction);
	// Update the coeficients
	superPair.a000 = coef[0];
	superPair.a001 = coef[1];
	superPair.a100 = coef[2];
	superPair.a101 = coef[3];
	superPair.a010 = coef[4];
	superPair.a011 = coef[5];
	superPair.a110 = coef[6];
	superPair.a111 = coef[7];

	// Add it to the list
	effCombiningList.push_front(superPair);

	return;
}

void NESuperCluster::participateIn(DissociationReaction& reaction,
		double *coef) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<NECluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Create a new SuperClusterDissociationPair
	SuperClusterDissociationPair superPair(
			&static_cast<NECluster&>(reaction.dissociating), &emittedCluster,
			&reaction);
	// Update the coeficients
	superPair.a00 = coef[0];
	superPair.a01 = coef[1];
	superPair.a10 = coef[2];
	superPair.a11 = coef[3];

	// Add it to the list
	effDissociatingList.push_front(superPair);

	return;
}

void NESuperCluster::emitFrom(DissociationReaction& reaction, double *coef) {

	// Create a new SuperClusterDissociationPair
	SuperClusterDissociationPair superPair(
			&static_cast<NECluster&>(reaction.first),
			&static_cast<NECluster&>(reaction.second), &reaction);
	// Update the coeficients
	superPair.a00 = coef[0];
	superPair.a01 = coef[1];
	superPair.a10 = coef[2];
	superPair.a11 = coef[3];

	// Add it to the list
	effEmissionList.push_front(superPair);

	return;
}

void NESuperCluster::updateFromNetwork() {

	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	// Aggregate the reacting pairs and combining reactants from the xeVector
	// Loop on the xeVector
	for (int i = 0; i < xeVector.size(); i++) {
		// Get the cluster size
		double size = xeVector[i]->getSize();
		// Get all vectors
		auto react = xeVector[i]->reactingPairs;
		auto combi = xeVector[i]->combiningReactants;
		auto disso = xeVector[i]->dissociatingPairs;
		auto emi = xeVector[i]->emissionPairs;

		// Set them in the super cluster map
		reactingMap[size] = react;
		combiningMap[size] = combi;
		dissociatingMap[size] = disso;
		emissionMap[size] = emi;
	}

	// Compute the dispersion
	computeDispersion();

	// Clear the xeVector
	xeVector.clear();

	return;
}

double NESuperCluster::getMoment() const {
	return l1;
}

double NESuperCluster::getTotalConcentration() const {
	// Initial declarations
	int index = 0;
	double distance = 0.0, conc = 0.0;

	// Loop on the xenon width
	for (int k = 0; k < sectionWidth; k++) {
		// Compute the xenon index
		index = (int) (numXe - (double) sectionWidth / 2.0) + k + 1;

		// Compute the distances
		distance = getDistance(index);

		// Add the concentration of each cluster in the group
		conc += getConcentration(distance);
	}

	return conc;
}

double NESuperCluster::getTotalXenonConcentration() const {
	// Initial declarations
	int index = 0;
	double distance = 0.0, conc = 0.0;

	// Loop on the xenon width
	for (int k = 0; k < sectionWidth; k++) {
		// Compute the xenon index
		index = (int) (numXe - (double) sectionWidth / 2.0) + k + 1;

		// Compute the distances
		distance = getDistance(index);

		// Add the concentration of each cluster in the group times the number of xenon atoms
		conc += getConcentration(distance) * (double) index;
	}

	return conc;
}

double NESuperCluster::getDistance(int xe) const {
	if (sectionWidth == 1)
		return 0.0;
	return 2.0 * (double) (xe - numXe) / ((double) sectionWidth - 1.0);
}

void NESuperCluster::computeDispersion() {
	// Local declarations
	int index = 0;
	// Initialize the dispersion sum
	double nXeSquare = 0.0;

	// Loop on the xenon width
	for (int k = 0; k < sectionWidth; k++) {
		// Compute the vacancy index
		index = (int) (numXe - (double) sectionWidth / 2.0) + k + 1;

		// Compute nSquare for the dispersion
		nXeSquare += (double) index * index;
	}

	// Compute the dispersion
	if (sectionWidth == 1)
		dispersion = 1.0;
	else {
		dispersion = 2.0
				* (nXeSquare
						- ((double) composition[toCompIdx(Species::Xe)]
								* ((double) composition[toCompIdx(Species::Xe)]
										/ (double) sectionWidth)))
				/ ((double) (sectionWidth * (sectionWidth - 1)));
	}

	return;
}

void NESuperCluster::optimizeReactions() {
	// Local declarations
	double factor = 0.0, distance = 0.0;
	NECluster *firstReactant, *secondReactant, *dissociatingCluster,
			*otherEmittedCluster, *firstCluster, *secondCluster;
	int index = 0;

	// Loop on the reacting map
	for (auto mapIt = reactingMap.begin(); mapIt != reactingMap.end();
			++mapIt) {
		// Get the pairs
		auto pairs = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = pairs.begin(); it != pairs.end();) {
			// Get the two reacting clusters
			firstReactant = (*it).first;
			secondReactant = (*it).second;

			// Create the corresponding production reaction
			std::unique_ptr<ProductionReaction> reaction(
					new ProductionReaction(*firstReactant, *secondReactant));
			// Add it to the network
			auto& prref = network.add(std::move(reaction));

			// Create a new SuperClusterProductionPair
			SuperClusterProductionPair superPair(firstReactant, secondReactant,
					&prref);

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != reactingMap.end();
					++mapItBis) {
				// Compute the xenon index
				index = mapItBis->first;
				factor = (double) (index - numXe) / dispersion;

				// Get the pairs
				auto pairsBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = pairsBis.begin(); itBis != pairsBis.end();) {
					// Get the two reacting clusters
					auto firstReactantBis = (*itBis).first;
					auto secondReactantBis = (*itBis).second;

					// Check if it is the same reaction
					if (firstReactantBis == firstReactant
							&& secondReactantBis == secondReactant) {
						superPair.a000 += 1.0;
						superPair.a001 += factor;
						superPair.a100 += (*itBis).firstDistance;
						superPair.a101 += (*itBis).firstDistance * factor;
						superPair.a010 += (*itBis).secondDistance;
						superPair.a011 += (*itBis).secondDistance * factor;
						superPair.a110 += (*itBis).firstDistance
								* (*itBis).secondDistance;
						superPair.a111 += (*itBis).firstDistance
								* (*itBis).secondDistance * factor;

						// Do not delete the element if it is the original one
						if (itBis == it) {
							++itBis;
							continue;
						}

						// Remove the reaction from the vector
						itBis = pairsBis.erase(itBis);
					}
					// Go to the next element
					else
						++itBis;
				}

				// Give back the pairs
				mapItBis->second = pairsBis;
			}

			// Add the super pair
			effReactingList.push_front(superPair);

			// Remove the reaction from the vector
			it = pairs.erase(it);
		}
	}

	// Loop on the combining map
	for (auto mapIt = combiningMap.begin(); mapIt != combiningMap.end();
			++mapIt) {
		// Get the pairs
		auto clusters = mapIt->second;

		// Loop over all the reacting pairs
		for (auto it = clusters.begin(); it != clusters.end();) {
			// Get the combining cluster
			NECluster& combiningReactant = *((*it).combining);

			// Create the corresponding production reaction
			std::unique_ptr<ProductionReaction> reaction(
					new ProductionReaction(*this, combiningReactant));
			// Add it to the network
			auto& prref = network.add(std::move(reaction));

			// Create a new SuperClusterProductionPair with NULL as the second cluster because
			// we do not need it
			SuperClusterProductionPair superPair(&combiningReactant, nullptr,
					&prref);

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != combiningMap.end();
					++mapItBis) {
				// Compute the xenon index
				index = mapItBis->first;
				distance = getDistance(index);
				factor = (double) (index - numXe) / dispersion;

				// Access the combining pairs
				auto& clustersBis = mapItBis->second;

				// Set the total number of reactants that produce to 
				// form this one.
				// May involve removing items from the vector of
				// combining cluster objects (clusersBis).
				// To avoid invalidating iterators into this vector,
				// we use the idiom of noting which need to get deleted,
				// a std::remove_if to move the doomed ones to the end
				// of the vector, and then erase to remove them.
				std::set<CombiningCluster*> doomedCombining;
				for (auto itBis = clustersBis.begin();
						itBis != clustersBis.end(); ++itBis) {

					// Access the current combining cluster
					NECluster& combiningReactantBis = *((*itBis).combining);

					// Check if it is the same reaction
					if (&combiningReactantBis == &combiningReactant) {
						superPair.a000 += 1.0;
						superPair.a001 += factor;
						superPair.a010 += (*itBis).distance;
						superPair.a011 += (*itBis).distance * factor;
						superPair.a100 += distance;
						superPair.a101 += distance * factor;
						superPair.a110 += (*itBis).distance * distance;
						superPair.a111 += (*itBis).distance * distance * factor;

						// Determine if we need to delete this item.
						// Do not delete if it is the original one.
						if (itBis != it) {
							// It is not the original one, so indicate it
							// needs to be removed.
							// NB: The expression with &(*iter) seems odd -
							// why not just use itBis?  We want the
							// address of the object, and itBis is an
							// iterator, not the address of the object itself.
							// So we dereference the iterator to access
							// the object, then take its address.
							doomedCombining.emplace(&(*itBis));
						}
					}
				}

				// Now that we know which combining clusters to delete,
				// Move them all to the end of the vector.
				auto firstToRemoveIter = std::remove_if(clustersBis.begin(),
						clustersBis.end(),
						[&doomedCombining](CombiningCluster& currCombining) {
							// See if currCombiningPair is in our set
							// of doomed items.
							auto diter = doomedCombining.find(&currCombining);
							return (diter != doomedCombining.end());
						});
				// Now that the doomed are all moved to be contiguous
				// and at the end of the vector, erase them.
				clustersBis.erase(firstToRemoveIter, clustersBis.end());
			}

			// Add the super pair
			effCombiningList.push_front(superPair);

			// Remove the reaction from the vector
			it = clusters.erase(it);
		}
	}

	// Loop on the dissociating map
	for (auto mapIt = dissociatingMap.begin(); mapIt != dissociatingMap.end();
			++mapIt) {
		// Get the pairs
		auto pairs = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = pairs.begin(); it != pairs.end();) {
			// Get the two reacting clusters
			dissociatingCluster = (*it).first;
			otherEmittedCluster = (*it).second;

			// Create a dissociation reaction
			std::unique_ptr<DissociationReaction> reaction(
					new DissociationReaction(*dissociatingCluster, *this,
							*otherEmittedCluster));
			// Add it to the network
			auto& drref = network.add(std::move(reaction));

			// Create a new SuperClusterDissociationPair
			SuperClusterDissociationPair superPair(dissociatingCluster,
					otherEmittedCluster, &drref);

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != dissociatingMap.end();
					++mapItBis) {
				// Compute the xenon index
				index = mapItBis->first;
				factor = (double) (index - numXe) / dispersion;

				// Get the pairs
				auto pairsBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = pairsBis.begin(); itBis != pairsBis.end();) {
					// Get the two reacting clusters
					auto dissociatingClusterBis = (*itBis).first;
					auto otherEmittedClusterBis = (*itBis).second;

					// Check if it is the same reaction
					if (dissociatingClusterBis == dissociatingCluster
							&& otherEmittedClusterBis == otherEmittedCluster) {
						superPair.a00 += 1.0;
						superPair.a01 += factor;
						superPair.a10 += (*itBis).firstDistance;
						superPair.a11 += (*itBis).firstDistance * factor;

						// Do not delete the element if it is the original one
						if (itBis == it) {
							++itBis;
							continue;
						}

						// Remove the reaction from the vector
						itBis = pairsBis.erase(itBis);
					}
					// Go to the next element
					else
						++itBis;
				}

				// Give back the pairs
				mapItBis->second = pairsBis;
			}

			// Add the super pair
			effDissociatingList.push_front(superPair);

			// Remove the reaction from the vector
			it = pairs.erase(it);
		}
	}

	// Loop on the emission map
	for (auto mapIt = emissionMap.begin(); mapIt != emissionMap.end();
			++mapIt) {
		// Get the pairs
		auto pairs = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = pairs.begin(); it != pairs.end();) {
			// Get the two reacting clusters
			firstCluster = (*it).first;
			secondCluster = (*it).second;

			// Create a dissociation reaction
			std::unique_ptr<DissociationReaction> reaction(
					new DissociationReaction(*this, *firstCluster,
							*secondCluster));
			// Add it to the network
			auto& drref = network.add(std::move(reaction));

			// Create a new SuperClusterDissociationPair
			SuperClusterDissociationPair superPair(firstCluster, secondCluster,
					&drref);

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != emissionMap.end();
					++mapItBis) {
				// Compute the xenon index
				index = mapItBis->first;
				distance = getDistance(index);
				factor = (double) (index - numXe) / dispersion;

				// Get the pairs
				auto pairsBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = pairsBis.begin(); itBis != pairsBis.end();) {
					// Get the two reacting clusters
					auto firstClusterBis = (*itBis).first;
					auto secondClusterBis = (*itBis).second;

					// Check if it is the same reaction
					if (firstClusterBis == firstCluster
							&& secondClusterBis == secondCluster) {
						superPair.a00 += 1.0;
						superPair.a01 += factor;
						superPair.a10 += distance;
						superPair.a11 += distance * factor;

						// Do not delete the element if it is the original one
						if (itBis == it) {
							++itBis;
							continue;
						}

						// Remove the reaction from the vector
						itBis = pairsBis.erase(itBis);
					}
					// Go to the next element
					else
						++itBis;
				}

				// Give back the pairs
				mapItBis->second = pairsBis;
			}

			// Add the super pair
			effEmissionList.push_front(superPair);

			// Remove the reaction from the vector
			it = pairs.erase(it);
		}
	}

	// Clear the maps because they won't be used anymore
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();
	reactingMap.clear();
	combiningMap.clear();
	dissociatingMap.clear();
	emissionMap.clear();
	return;
}

void NESuperCluster::resetConnectivities() {
	// Clear both sets
	reactionConnectivitySet.clear();
	dissociationConnectivitySet.clear();

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);
	setDissociationConnectivity(id);
	setReactionConnectivity(momId[0]);
	setDissociationConnectivity(momId[0]);

	// Loop over all the reacting pairs
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
		// The cluster is connecting to both clusters in the pair
		setReactionConnectivity((*it).first->getId());
		setReactionConnectivity((*it).first->getMomentId());
		setReactionConnectivity((*it).second->getId());
		setReactionConnectivity((*it).second->getMomentId());
	}

	// Loop over all the combining pairs
	for (auto it = effCombiningList.begin(); it != effCombiningList.end();
			++it) {
		// The cluster is connecting to the combining cluster
		setReactionConnectivity((*it).first->getId());
		setReactionConnectivity((*it).first->getMomentId());
	}

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end();
			++it) {
		// The cluster is connecting to the combining cluster
		setDissociationConnectivity((*it).first->getId());
		setDissociationConnectivity((*it).first->getMomentId());
	}

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	// Initialize the partial vector for the moment
	int dof = network.getDOF();
	momentPartials.resize(dof, 0.0);

	return;
}

double NESuperCluster::getTotalFlux(int i) {
	// Initialize the moment flux
	momentFlux = 0.0;

	// Get the fluxes
	double prodFlux = getProductionFlux(i);
	double dissFlux = getDissociationFlux(i);
	double combFlux = getCombinationFlux(i);
	double emissFlux = getEmissionFlux(i);

	return prodFlux - combFlux + dissFlux - emissFlux;
}

double NESuperCluster::getDissociationFlux(int xi) {
	// Initial declarations
	double flux = 0.0, value = 0.0;
	NECluster *dissociatingCluster = nullptr;

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end();
			++it) {
		// Get the dissociating clusters
		dissociatingCluster = (*it).first;
		double l0A = dissociatingCluster->getConcentration(0.0);
		double l1A = dissociatingCluster->getMoment();
		// Update the flux
		value = (*it).reaction.kConstant[xi] / (double) nTot;
		flux += value * ((*it).a00 * l0A + (*it).a10 * l1A);
		// Compute the moment fluxes
		momentFlux += value * ((*it).a01 * l0A + (*it).a11 * l1A);
	}

	// Return the flux
	return flux;
}

double NESuperCluster::getEmissionFlux(int xi) {
	// Initial declarations
	double flux = 0.0, value = 0.0;

	// Loop over all the emission pairs
	for (auto it = effEmissionList.begin(); it != effEmissionList.end(); ++it) {
		// Update the flux
		value = (*it).reaction.kConstant[xi] / (double) nTot;
		flux += value * ((*it).a00 * l0 + (*it).a10 * l1);
		// Compute the moment fluxes
		momentFlux -= value * ((*it).a01 * l0 + (*it).a11 * l1);
	}

	return flux;
}

double NESuperCluster::getProductionFlux(int xi) {
	// Local declarations
	double flux = 0.0, value = 0.0;
	NECluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Loop over all the reacting pairs
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
		// Get the two reacting clusters
		firstReactant = (*it).first;
		secondReactant = (*it).second;
		double l0A = firstReactant->getConcentration();
		double l0B = secondReactant->getConcentration();
		double l1A = firstReactant->getMoment();
		double l1B = secondReactant->getMoment();
		// Update the flux
		value = (*it).reaction.kConstant[xi] / (double) nTot;
		flux += value
				* ((*it).a000 * l0A * l0B + (*it).a010 * l0A * l1B
						+ (*it).a100 * l1A * l0B + (*it).a110 * l1A);
		// Compute the moment flux
		momentFlux += value
				* ((*it).a001 * l0A * l0B + (*it).a011 * l0A * l1B
						+ (*it).a101 * l1A * l0B + (*it).a111 * l1A);
	}

	// Return the production flux
	return flux;
}

double NESuperCluster::getCombinationFlux(int xi) {
	// Local declarations
	double flux = 0.0, value = 0.0;
	NECluster *combiningCluster = nullptr;

	// Loop over all the combining clusters
	for (auto it = effCombiningList.begin(); it != effCombiningList.end();
			++it) {
		// Get the two reacting clusters
		combiningCluster = (*it).first;
		double l0A = combiningCluster->getConcentration();
		double l1A = combiningCluster->getMoment();
		// Update the flux
		value = (*it).reaction.kConstant[xi] / (double) nTot;
		flux += value
				* ((*it).a000 * l0A * l0 + (*it).a100 * l0A * l1
						+ (*it).a010 * l1A * l0 + (*it).a110 * l1A * l1);
		// Compute the moment flux
		momentFlux -= value
				* ((*it).a001 * l0A * l0 + (*it).a101 * l0A * l1
						+ (*it).a011 * l1A * l0 + (*it).a111 * l1A * l1);
	}

	return flux;
}

void NESuperCluster::getPartialDerivatives(std::vector<double> & partials,
		int i) const {
	// Reinitialize the moment partial derivatives vector
	std::fill(momentPartials.begin(), momentPartials.end(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials, i);
	getCombinationPartialDerivatives(partials, i);
	getDissociationPartialDerivatives(partials, i);
	getEmissionPartialDerivatives(partials, i);

	return;
}

void NESuperCluster::getProductionPartialDerivatives(
		std::vector<double> & partials, int xi) const {
	// Initial declarations
	double value = 0.0;
	int index = 0;
	NECluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A

	// Loop over all the reacting pairs
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
		// Get the two reacting clusters
		firstReactant = (*it).first;
		secondReactant = (*it).second;
		double l0A = firstReactant->getConcentration();
		double l0B = secondReactant->getConcentration();
		double l1A = firstReactant->getMoment();
		double l1B = secondReactant->getMoment();

		// Compute the contribution from the first part of the reacting pair
		value = (*it).reaction.kConstant[xi] / (double) nTot;
		index = firstReactant->getId() - 1;
		partials[index] += value * ((*it).a000 * l0B + (*it).a010 * l1B);
		momentPartials[index] += value * ((*it).a001 * l0B + (*it).a011 * l1B);
		index = firstReactant->getMomentId() - 1;
		partials[index] += value * ((*it).a100 * l0B + (*it).a110 * l1B);
		momentPartials[index] += value * ((*it).a101 * l0B + (*it).a111 * l1B);
		// Compute the contribution from the second part of the reacting pair
		index = secondReactant->getId() - 1;
		partials[index] += value * ((*it).a000 * l0A + (*it).a100 * l1A);
		momentPartials[index] += value * ((*it).a001 * l0A + (*it).a101 * l1A);
		index = secondReactant->getMomentId() - 1;
		partials[index] += value * ((*it).a010 * l0A + (*it).a110 * l1A);
		momentPartials[index] += value * ((*it).a011 * l0A + (*it).a111 * l1A);
	}

	return;
}

void NESuperCluster::getCombinationPartialDerivatives(
		std::vector<double> & partials, int xi) const {
	// Initial declarations
	int index = 0;
	NECluster *cluster = nullptr;
	double value = 0.0;

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A

	// Loop over all the combining clusters
	for (auto it = effCombiningList.begin(); it != effCombiningList.end();
			++it) {
		// Get the two reacting clusters
		cluster = (*it).first;
		double l0A = cluster->getConcentration();
		double l1A = cluster->getMoment();

		// Compute the contribution from the combining cluster
		value = (*it).reaction.kConstant[xi] / (double) nTot;
		index = cluster->getId() - 1;
		partials[index] -= value * ((*it).a000 * l0 + (*it).a100 * l1);
		momentPartials[index] -= value * ((*it).a001 * l0 + (*it).a101 * l1);
		index = cluster->getMomentId() - 1;
		partials[index] -= value * ((*it).a010 * l0 + (*it).a110 * l1);
		momentPartials[index] -= value * ((*it).a011 * l0 + (*it).a111 * l1);
		// Compute the contribution from this cluster
		index = id - 1;
		partials[index] -= value * ((*it).a000 * l0A + (*it).a010 * l1A);
		momentPartials[index] -= value * ((*it).a001 * l0A + (*it).a011 * l1A);
		index = momId[0] - 1;
		partials[index] -= value * ((*it).a100 * l0A + (*it).a110 * l1A);
		momentPartials[index] -= value * ((*it).a101 * l0A + (*it).a111 * l1A);
	}

	return;
}

void NESuperCluster::getDissociationPartialDerivatives(
		std::vector<double> & partials, int xi) const {
	// Initial declarations
	int index = 0;
	NECluster *cluster = nullptr;
	double value = 0.0;

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end();
			++it) {
		// Get the dissociating clusters
		cluster = (*it).first;

		// Compute the contribution from the dissociating cluster
		value = (*it).reaction.kConstant[xi] / (double) nTot;
		index = cluster->getId() - 1;
		partials[index] += value * ((*it).a00);
		momentPartials[index] += value * ((*it).a01);
		index = cluster->getMomentId() - 1;
		partials[index] += value * ((*it).a10);
		momentPartials[index] += value * ((*it).a11);
	}

	return;
}

void NESuperCluster::getEmissionPartialDerivatives(
		std::vector<double> & partials, int xi) const {
	// Initial declarations
	int index = 0;
	double value = 0.0;

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)

	// Loop over all the emission pairs
	for (auto it = effEmissionList.begin(); it != effEmissionList.end(); ++it) {
		// Compute the contribution from the dissociating cluster
		value = (*it).reaction.kConstant[xi] / (double) nTot;
		index = id - 1;
		partials[index] -= value * ((*it).a00);
		momentPartials[index] -= value * ((*it).a01);
		index = momId[0] - 1;
		partials[index] -= value * ((*it).a10);
		momentPartials[index] -= value * ((*it).a11);
	}

	return;
}

void NESuperCluster::getMomentPartialDerivatives(
		std::vector<double> & partials) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = momentPartials[i];
	}

	return;
}

std::vector<std::vector<double> > NESuperCluster::getProdVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the reacting pairs
	std::for_each(effReactingList.begin(), effReactingList.end(),
			[&toReturn](SuperClusterProductionPair const& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.second->getId() - 1);
				tempVec.push_back(currPair.a000);
				tempVec.push_back(currPair.a001);
				tempVec.push_back(currPair.a100);
				tempVec.push_back(currPair.a101);
				tempVec.push_back(currPair.a010);
				tempVec.push_back(currPair.a011);
				tempVec.push_back(currPair.a110);
				tempVec.push_back(currPair.a111);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > NESuperCluster::getCombVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the combining reactants
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[&toReturn](SuperClusterProductionPair const& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.a000);
				tempVec.push_back(currPair.a001);
				tempVec.push_back(currPair.a100);
				tempVec.push_back(currPair.a101);
				tempVec.push_back(currPair.a010);
				tempVec.push_back(currPair.a011);
				tempVec.push_back(currPair.a110);
				tempVec.push_back(currPair.a111);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > NESuperCluster::getDissoVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the dissociating pairs
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[&toReturn](SuperClusterDissociationPair const& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.second->getId() - 1);
				tempVec.push_back(currPair.a00);
				tempVec.push_back(currPair.a01);
				tempVec.push_back(currPair.a10);
				tempVec.push_back(currPair.a11);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > NESuperCluster::getEmitVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the emitting pairs
	std::for_each(effEmissionList.begin(), effEmissionList.end(),
			[&toReturn](SuperClusterDissociationPair const& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.second->getId() - 1);
				tempVec.push_back(currPair.a00);
				tempVec.push_back(currPair.a01);
				tempVec.push_back(currPair.a10);
				tempVec.push_back(currPair.a11);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}
