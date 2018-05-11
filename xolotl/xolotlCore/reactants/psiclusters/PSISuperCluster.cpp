// Includes
#include <iterator>
#include "PSISuperCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <MathUtils.h>
#include <xolotlPerf.h>

using namespace xolotlCore;

PSISuperCluster::PSISuperCluster(double num[4], int _nTot, int width[4],
		IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(_network, registry,
				buildName(num[0], num[1], num[2], num[3])), nTot(_nTot), l0(0.0) {
	// Loop on the axis
	for (int i = 0; i < 4; i++) {
		// Set the cluster size as the sum of
		// the number of Helium and Vacancies
		numAtom[i] = num[i];
		size += (int) num[i];
		// Set the width
		sectionWidth[i] = width[i];
	}

	// Update the composition map
	composition[toCompIdx(Species::He)] = (int) num[0];
	composition[toCompIdx(Species::D)] = (int) num[1];
	composition[toCompIdx(Species::T)] = (int) num[2];
	composition[toCompIdx(Species::V)] = (int) num[3];

	// Set the formation energy
	formationEnergy = 0.0; // It is set to 0.0 because we do not want the super clusters to undergo dissociation

	// Set the diffusion factor and the migration energy
	migrationEnergy = std::numeric_limits<double>::infinity();
	diffusionFactor = 0.0;

	// Set the typename appropriately
	type = ReactantType::PSISuper;

	return;
}

void PSISuperCluster::resultFrom(ProductionReaction& reaction, int a[4],
		int b[4]) {

	// Check if we already know about the reaction.
	auto rkey = std::make_pair(&(reaction.first), &(reaction.second));
	auto it = effReactingList.find(rkey);
	if (it == effReactingList.end()) {

		// We did not already know about this reaction.
		// Add info about production to our list.
		auto eret = effReactingList.emplace(std::piecewise_construct,
				std::forward_as_tuple(rkey),
				std::forward_as_tuple(reaction,
						static_cast<PSICluster&>(reaction.first),
						static_cast<PSICluster&>(reaction.second)));
		// Since we already checked and didn't know about the reaction,
		// we had better have added it with our emplace() call.
		assert(eret.second);
		it = eret.first;
	}
	assert(it != effReactingList.end());
	auto& prodPair = it->second;

	// NB: prodPair's reactants are same as reaction.
	// So use prodPair only from here on.
	// TODO any way to enforce this?

	// Update the coefficients
	double firstDistance[5] = { }, secondDistance[5] = { };
	firstDistance[0] = 1.0, secondDistance[0] = 1.0;
	if (prodPair.first.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(prodPair.first);
		for (int i = 1; i < psDim; i++) {
			firstDistance[i] = super.getDistance(b[indexList[i] - 1],
					indexList[i] - 1);
		}
	}
	if (prodPair.second.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(prodPair.second);
		for (int i = 1; i < psDim; i++) {
			secondDistance[i] = super.getDistance(b[indexList[i] - 1],
					indexList[i] - 1);
		}
	}
	double factor[5] = { };
	factor[0] = 1.0;
	for (int i = 1; i < psDim; i++) {
		factor[i] = getFactor(a[indexList[i] - 1], indexList[i] - 1);
	}
	// First is A, second is B, in A + B -> this
	for (int k = 0; k < psDim; k++) {
		for (int j = 0; j < psDim; j++) {
			for (int i = 0; i < psDim; i++) {
				prodPair.coefs[i][j][k] += firstDistance[i] * secondDistance[j]
						* factor[k];
			}
		}
	}

	return;
}

void PSISuperCluster::resultFrom(ProductionReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& prInfos) {

	// Check if we already know about the reaction.
	auto rkey = std::make_pair(&(reaction.first), &(reaction.second));
	auto it = effReactingList.find(rkey);
	if (it == effReactingList.end()) {

		// We did not already know about this reaction.
		// Add info about production to our list.
		auto eret = effReactingList.emplace(std::piecewise_construct,
				std::forward_as_tuple(rkey),
				std::forward_as_tuple(reaction,
						static_cast<PSICluster&>(reaction.first),
						static_cast<PSICluster&>(reaction.second)));
		// Since we already checked and didn't know about the reaction,
		// we had better have added it with our emplace() call.
		assert(eret.second);
		it = eret.first;
	}
	assert(it != effReactingList.end());
	auto& prodPair = it->second;

	// NB: prodPair's reactants are same as reaction.
	// So use prodPair only from here on.
	// TODO any way to enforce this?

	// Update the coefficients
	std::for_each(prInfos.begin(), prInfos.end(),
			[this,&prodPair](const PendingProductionReactionInfo& currPRI) {
				// Update the coefficients
				double firstDistance[5] = {}, secondDistance[5] = {};
				firstDistance[0] = 1.0, secondDistance[0] = 1.0;
				if (prodPair.first.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster const&>(prodPair.first);
					for (int i = 1; i < psDim; i++) {
						firstDistance[i] = super.getDistance(currPRI.b[indexList[i]-1], indexList[i]-1);
					}
				}
				if (prodPair.second.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster const&>(prodPair.second);
					for (int i = 1; i < psDim; i++) {
						secondDistance[i] = super.getDistance(currPRI.b[indexList[i]-1], indexList[i]-1);
					}
				}
				double factor[5] = {};
				factor[0] = 1.0;
				for (int i = 1; i < psDim; i++) {
					factor[i] = getFactor(currPRI.a[indexList[i]-1], indexList[i]-1);
				}
				// First is A, second is B, in A + B -> this
				for (int k = 0; k < psDim; k++) {
					for (int j = 0; j < psDim; j++) {
						for (int i = 0; i < psDim; i++) {
							prodPair.coefs[i][j][k] += firstDistance[i] * secondDistance[j]
							* factor[k];
						}
					}
				}
			});

	return;
}

void PSISuperCluster::participateIn(ProductionReaction& reaction, int a[4]) {

	setReactionConnectivity(id);
	// Look for the other cluster
	auto& otherCluster = static_cast<PSICluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if we already know about the reaction.
	auto rkey = &otherCluster;
	auto it = effCombiningList.find(rkey);
	if (it == effCombiningList.end()) {

		// We did not already know about the reaction.
		// Note that we combine with the other cluster in this reaction.
		auto eret = effCombiningList.emplace(std::piecewise_construct,
				std::forward_as_tuple(rkey),
				std::forward_as_tuple(reaction,
						static_cast<PSICluster&>(otherCluster)));
		// Since we already checked and didn't know about the reaction then,
		// we had better have added it with our emplace call.
		assert(eret.second);
		it = eret.first;
	}
	assert(it != effCombiningList.end());
	auto& combCluster = it->second;

	// Update the coefficients
	double distance[5] = { }, factor[5] = { };
	distance[0] = 1.0, factor[0] = 1.0;
	for (int i = 1; i < psDim; i++) {
		distance[i] = getDistance(a[indexList[i] - 1], indexList[i] - 1);
		factor[i] = getFactor(a[indexList[i] - 1], indexList[i] - 1);
	}

	// This is A, itBis is B, in A + B -> C
	for (int k = 0; k < psDim; k++) {
		for (int j = 0; j < psDim; j++) {
			combCluster.coefs[j][0][k] += distance[j] * factor[k];
		}
	}

	return;
}

void PSISuperCluster::participateIn(ProductionReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& pendingPRInfos) {

	setReactionConnectivity(id);
	// Look for the other cluster
	auto& otherCluster = static_cast<PSICluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if we already know about the reaction.
	auto rkey = &otherCluster;
	auto it = effCombiningList.find(rkey);
	if (it == effCombiningList.end()) {

		// We did not already know about the reaction.
		// Note that we combine with the other cluster in this reaction.
		auto eret = effCombiningList.emplace(std::piecewise_construct,
				std::forward_as_tuple(rkey),
				std::forward_as_tuple(reaction,
						static_cast<PSICluster&>(otherCluster)));
		// Since we already checked and didn't know about the reaction then,
		// we had better have added it with our emplace call.
		assert(eret.second);
		it = eret.first;
	}
	assert(it != effCombiningList.end());
	auto& combCluster = it->second;

	// Update the coefficients
	std::for_each(pendingPRInfos.begin(), pendingPRInfos.end(),
			[this,&combCluster,&otherCluster](const PendingProductionReactionInfo& currPRInfo) {
				// Update the coefficients
				double distance[5] = {}, factor[5] = {};
				distance[0] = 1.0, factor[0] = 1.0;
				for (int i = 1; i < psDim; i++) {
					distance[i] = getDistance(currPRInfo.b[indexList[i]-1], indexList[i]-1);
					factor[i] = getFactor(currPRInfo.b[indexList[i]-1], indexList[i]-1);
				}

				// This is A, itBis is B, in A + B -> C
				for (int k = 0; k < psDim; k++) {
					for (int j = 0; j < psDim; j++) {
						combCluster.coefs[j][0][k] += distance[j] * factor[k];
					}
				}
			});

	return;
}

void PSISuperCluster::participateIn(DissociationReaction& reaction, int a[4],
		int b[4]) {

	// Determine which is the other cluster.
	auto& emittedCluster = static_cast<PSICluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if we already know about the reaction.
	auto rkey = std::make_pair(&(reaction.dissociating), &emittedCluster);
	auto it = effDissociatingList.find(rkey);
	if (it == effDissociatingList.end()) {

		// We did not already know about it.

		// Add it to the network
		auto eret = effDissociatingList.emplace(std::piecewise_construct,
				std::forward_as_tuple(rkey),
				std::forward_as_tuple(reaction,
						static_cast<PSICluster&>(reaction.dissociating),
						static_cast<PSICluster&>(emittedCluster)));
		// Since we already checked and didn't know about the reaction then,
		// we had better have added it with our emplace() call.
		assert(eret.second);
		it = eret.first;
	}
	assert(it != effDissociatingList.end());
	auto& dissPair = it->second;

	// Update the coefficients
	double distance[5] = { }, factor[5] = { };
	distance[0] = 1.0, factor[0] = 1.0;
	if (reaction.dissociating.getType() == ReactantType::PSISuper) {
		auto const& super =
				static_cast<PSICluster const&>(reaction.dissociating);
		for (int i = 1; i < psDim; i++) {
			distance[i] = super.getDistance(a[indexList[i] - 1],
					indexList[i] - 1);
		}
	}
	for (int i = 1; i < psDim; i++) {
		factor[i] = getFactor(b[indexList[i] - 1], indexList[i] - 1);
	}
	// A is the dissociating cluster
	for (int j = 0; j < psDim; j++) {
		for (int i = 0; i < psDim; i++) {
			dissPair.coefs[i][j] += distance[i] * factor[j];
		}
	}

	return;
}

void PSISuperCluster::participateIn(DissociationReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& prInfos) {

	// Determine which is the other cluster.
	auto& emittedCluster = static_cast<PSICluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if we already know about the reaction.
	auto rkey = std::make_pair(&(reaction.dissociating), &emittedCluster);
	auto it = effDissociatingList.find(rkey);
	if (it == effDissociatingList.end()) {

		// We did not already know about it.

		// Add it to the network
		auto eret = effDissociatingList.emplace(std::piecewise_construct,
				std::forward_as_tuple(rkey),
				std::forward_as_tuple(reaction,
						static_cast<PSICluster&>(reaction.dissociating),
						static_cast<PSICluster&>(emittedCluster)));
		// Since we already checked and didn't know about the reaction then,
		// we had better have added it with our emplace() call.
		assert(eret.second);
		it = eret.first;
	}
	assert(it != effDissociatingList.end());
	auto& dissPair = it->second;

	// Update the coefficients
	std::for_each(prInfos.begin(), prInfos.end(),
			[this,&dissPair,&reaction](const PendingProductionReactionInfo& currPRI) {
				// Update the coefficients
				double distance[5] = {}, factor[5] = {};
				distance[0] = 1.0, factor[0] = 1.0;
				if (reaction.dissociating.getType() == ReactantType::PSISuper) {
					auto const& super =
					static_cast<PSICluster const&>(reaction.dissociating);
					for (int i = 1; i < psDim; i++) {
						distance[i] = super.getDistance(currPRI.a[indexList[i]-1], indexList[i]-1);
					}
				}
				for (int i = 1; i < psDim; i++) {
					factor[i] = getFactor(currPRI.b[indexList[i]-1], indexList[i]-1);
				}
				// A is the dissociating cluster
				for (int j = 0; j < psDim; j++) {
					for (int i = 0; i < psDim; i++) {
						dissPair.coefs[i][j] += distance[i] * factor[j];
					}
				}
			});

	return;
}

void PSISuperCluster::emitFrom(DissociationReaction& reaction, int a[4]) {

	// Check if we already know about the reaction.
	auto rkey = std::make_pair(&(reaction.first), &(reaction.second));
	auto it = effEmissionList.find(rkey);
	if (it == effEmissionList.end()) {

		// We did not already know about it.

		// Note that we emit from the two rectants according to the given
		// reaction.
		auto eret = effEmissionList.emplace(std::piecewise_construct,
				std::forward_as_tuple(rkey),
				std::forward_as_tuple(reaction,
						static_cast<PSICluster&>(reaction.first),
						static_cast<PSICluster&>(reaction.second)));
		// Since we already checked and didn't know about the reaction then,
		// we had better have added it with our emplace() call.
		assert(eret.second);
		it = eret.first;
	}
	assert(it != effEmissionList.end());
	auto& dissPair = it->second;

	// Update the coefficients
	double distance[5] = { }, factor[5] = { };
	distance[0] = 1.0, factor[0] = 1.0;
	for (int i = 1; i < psDim; i++) {
		distance[i] = getDistance(a[indexList[i] - 1], indexList[i] - 1);
		factor[i] = getFactor(a[indexList[i] - 1], indexList[i] - 1);
	}
	// A is the dissociating cluster
	for (int j = 0; j < psDim; j++) {
		for (int i = 0; i < psDim; i++) {
			dissPair.coefs[i][j] += distance[i] * factor[j];
		}
	}

	return;
}

void PSISuperCluster::emitFrom(DissociationReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& prInfos) {

	// Check if we already know about the reaction.
	auto rkey = std::make_pair(&(reaction.first), &(reaction.second));
	auto it = effEmissionList.find(rkey);
	if (it == effEmissionList.end()) {

		// We did not already know about it.

		// Note that we emit from the two rectants according to the given
		// reaction.
		auto eret = effEmissionList.emplace(std::piecewise_construct,
				std::forward_as_tuple(rkey),
				std::forward_as_tuple(reaction,
						static_cast<PSICluster&>(reaction.first),
						static_cast<PSICluster&>(reaction.second)));
		// Since we already checked and didn't know about the reaction then,
		// we had better have added it with our emplace() call.
		assert(eret.second);
		it = eret.first;
	}
	assert(it != effEmissionList.end());
	auto& dissPair = it->second;

	// Update the coefficients
	std::for_each(prInfos.begin(), prInfos.end(),
			[this,&dissPair](const PendingProductionReactionInfo& currPRI) {
				// Update the coefficients
				double distance[5] = {}, factor[5] = {};
				distance[0] = 1.0, factor[0] = 1.0;
				for (int i = 1; i < psDim; i++) {
					distance[i] = getDistance(currPRI.a[indexList[i]-1], indexList[i]-1);
					factor[i] = getFactor(currPRI.a[indexList[i]-1], indexList[i]-1);
				}
				// A is the dissociating cluster
				for (int j = 0; j < psDim; j++) {
					for (int i = 0; i < psDim; i++) {
						dissPair.coefs[i][j] += distance[i] * factor[j];
					}
				}
			});

	return;
}

void PSISuperCluster::setHeVVector(
		std::set<std::tuple<int, int, int, int> > vec) {
	// Copy the list of coordinates
	heVList = vec;

	// Initialize the dispersion sum
	double nSquare[4] = { };
	// Update the network map, compute the radius and dispersions
	for (auto const& pair : heVList) {
		double rad = (sqrt(3.0) / 4.0) * xolotlCore::tungstenLatticeConstant
				+ pow(
						(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0)
								* std::get<3>(pair)) / (8.0 * xolotlCore::pi),
						(1.0 / 3.0))
				- pow(
						(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0))
								/ (8.0 * xolotlCore::pi), (1.0 / 3.0));
		reactionRadius += rad / (double) nTot;

		// Compute nSquare for the dispersion
		nSquare[0] += (double) (std::get<0>(pair) * std::get<0>(pair));
		nSquare[1] += (double) (std::get<1>(pair) * std::get<1>(pair));
		nSquare[2] += (double) (std::get<2>(pair) * std::get<2>(pair));
		nSquare[3] += (double) (std::get<3>(pair) * std::get<3>(pair));
	}

	// Loop on the different type of clusters in grouping
	for (int i = 0; i < 4; i++) {
		// Compute the dispersions
		if (sectionWidth[i] == 1)
			dispersion[i] = 1.0;
		else
			dispersion[i] = 2.0
					* (nSquare[i] - (numAtom[i] * (double) nTot * numAtom[i]))
					/ ((double) (nTot * (sectionWidth[i] - 1)));

		// Set the boundaries
		bounds[i] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>((numAtom[i]
						- (double) sectionWidth[i] / 2.0) + 1),
				static_cast<IReactant::SizeType>((numAtom[i]
						- (double) sectionWidth[i] / 2.0) + sectionWidth[i])
						+ 1);
	}

	return;
}

double PSISuperCluster::getTotalConcentration() const {
	// Initial declarations
	double heDistance = 0.0, dDistance = 0.0, tDistance = 0.0, vDistance = 0.0,
			conc = 0.0;

	// Loop on the indices
	for (auto const& pair : heVList) {
		// Compute the distances
		heDistance = getDistance(std::get<0>(pair), 0);
		dDistance = getDistance(std::get<1>(pair), 1);
		tDistance = getDistance(std::get<2>(pair), 2);
		vDistance = getDistance(std::get<3>(pair), 3);

		// Add the concentration of each cluster in the group times its number of helium
		conc += getConcentration(heDistance, dDistance, tDistance, vDistance);
	}

	return conc;
}

double PSISuperCluster::getTotalAtomConcentration(int axis) const {
	// Initial declarations
	double heDistance = 0.0, dDistance = 0.0, tDistance = 0.0, vDistance = 0.0,
			conc = 0.0;

	if (axis == 0) {
		// Loop on the indices
		for (auto const& pair : heVList) {
			// Compute the distances
			heDistance = getDistance(std::get<0>(pair), 0);
			dDistance = getDistance(std::get<1>(pair), 1);
			tDistance = getDistance(std::get<2>(pair), 2);
			vDistance = getDistance(std::get<3>(pair), 3);

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, dDistance, tDistance,
					vDistance) * (double) std::get<0>(pair);
		}
	} else if (axis == 1) {
		// Loop on the indices
		for (auto const& pair : heVList) {
			// Compute the distances
			heDistance = getDistance(std::get<0>(pair), 0);
			dDistance = getDistance(std::get<1>(pair), 1);
			tDistance = getDistance(std::get<2>(pair), 2);
			vDistance = getDistance(std::get<3>(pair), 3);

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, dDistance, tDistance,
					vDistance) * (double) std::get<1>(pair);
		}
	} else if (axis == 2) {
		// Loop on the indices
		for (auto const& pair : heVList) {
			// Compute the distances
			heDistance = getDistance(std::get<0>(pair), 0);
			dDistance = getDistance(std::get<1>(pair), 1);
			tDistance = getDistance(std::get<2>(pair), 2);
			vDistance = getDistance(std::get<3>(pair), 3);

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, dDistance, tDistance,
					vDistance) * (double) std::get<2>(pair);
		}
	}

	return conc;
}

double PSISuperCluster::getTotalVacancyConcentration() const {
	// Initial declarations
	double heDistance = 0.0, dDistance = 0.0, tDistance = 0.0, vDistance = 0.0,
			conc = 0.0;

	// Loop on the indices
	for (auto const& pair : heVList) {
		// Compute the distances
		heDistance = getDistance(std::get<0>(pair), 0);
		dDistance = getDistance(std::get<1>(pair), 1);
		tDistance = getDistance(std::get<2>(pair), 2);
		vDistance = getDistance(std::get<3>(pair), 3);

		// Add the concentration of each cluster in the group times its number of vacancy
		conc += getConcentration(heDistance, dDistance, tDistance, vDistance)
				* (double) std::get<3>(pair);
	}

	return conc;
}

double PSISuperCluster::getIntegratedVConcentration(int v) const {
	// Initial declarations
	double heDistance = 0.0, dDistance = 0.0, tDistance = 0.0, vDistance = 0.0,
			conc = 0.0;

	// Loop on the indices
	for (auto const& pair : heVList) {
		// Skip the wrong V size
		if (std::get<3>(pair) != v)
			continue;

		// Compute the distances
		heDistance = getDistance(std::get<0>(pair), 0);
		dDistance = getDistance(std::get<1>(pair), 1);
		tDistance = getDistance(std::get<2>(pair), 2);
		vDistance = getDistance(std::get<3>(pair), 3);

		// Add the concentration of each cluster in the group
		conc += getConcentration(heDistance, dDistance, tDistance, vDistance);
	}

	return conc;
}

void PSISuperCluster::resetConnectivities() {
	// Clear both sets
	reactionConnectivitySet.clear();
	dissociationConnectivitySet.clear();

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);
	setDissociationConnectivity(id);
	for (int i = 0; i < 4; i++) {
		setReactionConnectivity(momId[i]);
		setDissociationConnectivity(momId[i]);
	}
	// Visit all the reacting pairs
	std::for_each(effReactingList.begin(), effReactingList.end(),
			[this](ProductionPairMap::value_type const& currMapItem) {
				// The cluster is connecting to both clusters in the pair
				auto const& currPair = currMapItem.second;
				setReactionConnectivity(currPair.first.getId());
				setReactionConnectivity(currPair.second.getId());
				for (int i = 1; i < psDim; i++) {
					setReactionConnectivity(currPair.first.getMomentId(indexList[i]-1));
					setReactionConnectivity(currPair.second.getMomentId(indexList[i]-1));
				}
			});

	// Visit all the combining pairs
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[this](CombiningClusterMap::value_type const& currMapItem) {
				// The cluster is connecting to the combining cluster
				auto const& currComb = currMapItem.second;
				setReactionConnectivity(currComb.first.getId());
				for (int i = 1; i < psDim; i++) {
					setReactionConnectivity(currComb.first.getMomentId(indexList[i]-1));
				}
			});

	// Loop over all the dissociating pairs
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this](DissociationPairMap::value_type const& currMapItem) {
				// The cluster is connecting to the combining cluster
				auto const& currPair = currMapItem.second;
				setDissociationConnectivity(currPair.first.getId());
				for (int i = 1; i < psDim; i++) {
					setDissociationConnectivity(currPair.first.getMomentId(indexList[i]-1));
				}
			});

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	return;
}

double PSISuperCluster::getDissociationFlux() {
	// Initial declarations
	double flux = 0.0;

	// Sum over all the dissociating pairs
	// TODO consider using std::accumulate.  May also want to change side
	// effect of updating member variables heMomentFlux and
	// vMomentFlux here.
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this,&flux](DissociationPairMap::value_type const& currMapItem) {
				auto const& currPair = currMapItem.second;

				// Get the dissociating clusters
				auto const& dissociatingCluster = currPair.first;
				double lA[5] = {};
				lA[0] = dissociatingCluster.getConcentration();
				for (int i = 1; i < psDim; i++) {
					lA[i] = dissociatingCluster.getMoment(indexList[i]-1);
				}

				double sum[5] = {};
				for (int j = 0; j < psDim; j++) {
					for (int i = 0; i < psDim; i++) {
						sum[j] += currPair.coefs[i][j] * lA[i];
					}
				}
				// Update the flux
				auto value = currPair.kConstant / (double) nTot;
				flux += value * sum[0];
				// Compute the moment fluxes
				for (int i = 1; i < psDim; i++) {
					momentFlux[indexList[i]-1] += value * sum[i];
				}
			});

	// Return the flux
	return flux;
}

double PSISuperCluster::getEmissionFlux() {
	// Initial declarations
	double flux = 0.0;

	// Loop over all the emission pairs
	// TODO consider using std::accumulate.  May also want to change side
	// effect of updating member variables heMomentFlux and
	// vMomentFlux here.
	std::for_each(effEmissionList.begin(), effEmissionList.end(),
			[this,&flux](DissociationPairMap::value_type const& currMapItem) {
				auto const& currPair = currMapItem.second;
				double lA[5] = {};
				lA[0] = l0;
				for (int i = 1; i < psDim; i++) {
					lA[i] = l1[indexList[i]-1];
				}

				double sum[5] = {};
				for (int j = 0; j < psDim; j++) {
					for (int i = 0; i < psDim; i++) {
						sum[j] += currPair.coefs[i][j] * lA[i];
					}
				}
				// Update the flux
				auto value = currPair.kConstant / (double) nTot;
				flux += value * sum[0];
				// Compute the moment fluxes
				for (int i = 1; i < psDim; i++) {
					momentFlux[indexList[i]-1] -= value * sum[i];
				}
			});

	return flux;
}

double PSISuperCluster::getProductionFlux() {
	// Local declarations
	double flux = 0.0;

	// Sum over all the reacting pairs
	// TODO consider using std::accumulate.  May also want to change side
	// effect of updating member variables heMomentFlux and
	// vMomentFlux here.
	std::for_each(effReactingList.begin(), effReactingList.end(),
			[this,&flux](ProductionPairMap::value_type const& currMapItem) {

				auto const& currPair = currMapItem.second;

				// Get the two reacting clusters
				auto const& firstReactant = currPair.first;
				auto const& secondReactant = currPair.second;
				double lA[5] = {}, lB[5] = {};
				lA[0] = firstReactant.getConcentration();
				lB[0] = secondReactant.getConcentration();
				for (int i = 1; i < psDim; i++) {
					lA[i] = firstReactant.getMoment(indexList[i]-1);
					lB[i] = secondReactant.getMoment(indexList[i]-1);
				}

				double sum[5] = {};
				for (int k = 0; k < psDim; k++) {
					for (int j = 0; j < psDim; j++) {
						for (int i = 0; i < psDim; i++) {
							sum[k] += currPair.coefs[j][i][k] * lA[j] * lB[i];
						}
					}
				}

				// Update the flux
				auto value = currPair.kConstant / (double) nTot;
				flux += value * sum[0];
				// Compute the moment fluxes
				for (int i = 1; i < psDim; i++) {
					momentFlux[indexList[i]-1] += value * sum[i];
				}
			});

	// Return the production flux
	return flux;
}

double PSISuperCluster::getCombinationFlux() {
	// Local declarations
	double flux = 0.0;

	// Sum over all the combining clusters
	// TODO consider using std::accumulate.  May also want to change side
	// effect of updating member variables heMomentFlux and
	// vMomentFlux here.
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[this,&flux](CombiningClusterMap::value_type const& currMapItem) {
				// Get the combining cluster
				auto const& currComb = currMapItem.second;
				auto const& combiningCluster = currComb.first;
				double lA[5] = {}, lB[5] = {};
				lA[0] = l0;
				lB[0] = combiningCluster.getConcentration();
				for (int i = 1; i < psDim; i++) {
					lA[i] = l1[indexList[i]-1];
					lB[i] = combiningCluster.getMoment(indexList[i]-1);
				}

				double sum[5] = {};
				for (int k = 0; k < psDim; k++) {
					for (int j = 0; j < psDim; j++) {
						for (int i = 0; i < psDim; i++) {
							sum[k] += currComb.coefs[i][j][k] * lA[i] * lB[j];
						}
					}
				}
				// Update the flux
				auto value = currComb.kConstant / (double) nTot;
				flux += value * sum[0];
				// Compute the moment fluxes
				for (int i = 1; i < psDim; i++) {
					momentFlux[indexList[i]-1] -= value * sum[i];
				}
			});

	return flux;
}

void PSISuperCluster::computePartialDerivatives(double* partials[5],
		const ReactionNetwork::PartialsIdxMap partialsIdxMap[5]) const {

	// Get the partial derivatives for each reaction type
	computeProductionPartialDerivatives(partials, partialsIdxMap);
	computeCombinationPartialDerivatives(partials, partialsIdxMap);
	computeDissociationPartialDerivatives(partials, partialsIdxMap);
	computeEmissionPartialDerivatives(partials, partialsIdxMap);

	return;
}

void PSISuperCluster::computeProductionPartialDerivatives(double* partials[5],
		const ReactionNetwork::PartialsIdxMap partialsIdxMap[5]) const {

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A

	// Loop over all the reacting pairs
	std::for_each(effReactingList.begin(), effReactingList.end(),
			[this,
			&partials, &partialsIdxMap](ProductionPairMap::value_type const& currMapItem) {

				auto const& currPair = currMapItem.second;
				// Get the two reacting clusters
				auto const& firstReactant = currPair.first;
				auto const& secondReactant = currPair.second;
				double lA[5] = {}, lB[5] = {};
				lA[0] = firstReactant.getConcentration();
				lB[0] = secondReactant.getConcentration();
				for (int i = 1; i < psDim; i++) {
					lA[i] = firstReactant.getMoment(indexList[i]-1);
					lB[i] = secondReactant.getMoment(indexList[i]-1);
				}

				double sum[5][5][2] = {};
				for (int k = 0; k < psDim; k++) {
					for (int j = 0; j < psDim; j++) {
						for (int i = 0; i < psDim; i++) {
							sum[k][j][0] += currPair.coefs[j][i][k] * lB[i];
							sum[k][j][1] += currPair.coefs[i][j][k] * lA[i];
						}
					}
				}

				// Compute the contribution from the first and second part of the reacting pair
				auto value = currPair.kConstant / (double) nTot;
				for (int j = 0; j < psDim; j++) {
					int indexA = 0, indexB = 0;
					if (j == 0) {
						indexA = firstReactant.getId() - 1;
						indexB = secondReactant.getId() - 1;
					}
					else {
						indexA = firstReactant.getMomentId(indexList[j]-1) - 1;
						indexB = secondReactant.getMomentId(indexList[j]-1) - 1;
					}
					auto partialsIdxA = partialsIdxMap[j].at(indexA);
					auto partialsIdxB = partialsIdxMap[j].at(indexB);
					for (int i = 0; i < psDim; i++) {
						partials[i][partialsIdxA] += value * sum[i][j][0];
						partials[i][partialsIdxB] += value * sum[i][j][1];
					}
				}
			});

	return;
}

void PSISuperCluster::computeCombinationPartialDerivatives(double* partials[5],
		const ReactionNetwork::PartialsIdxMap partialsIdxMap[5]) const {

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A

	// Visit all the combining clusters
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[this,
			&partials, &partialsIdxMap](CombiningClusterMap::value_type const& currMapItem) {
				// Get the combining clusters
				auto const& currComb = currMapItem.second;
				auto const& cluster = currComb.first;
				double lA[5] = {}, lB[5] = {};
				lA[0] = l0;
				lB[0] = cluster.getConcentration();
				for (int i = 1; i < psDim; i++) {
					lA[i] = l1[indexList[i]-1];
					lB[i] = cluster.getMoment(indexList[i]-1);
				}

				double sum[5][5][2] = {};
				for (int k = 0; k < psDim; k++) {
					for (int j = 0; j < psDim; j++) {
						for (int i = 0; i < psDim; i++) {
							sum[k][j][0] += currComb.coefs[i][j][k] * lA[i];
							sum[k][j][1] += currComb.coefs[j][i][k] * lB[i];
						}
					}
				}

				// Compute the contribution from the both clusters
				auto value = currComb.kConstant / (double) nTot;
				for (int j = 0; j < psDim; j++) {
					int indexA = 0, indexB = 0;
					if (j == 0) {
						indexA = cluster.getId() - 1;
						indexB = id - 1;
					}
					else {
						indexA = cluster.getMomentId(indexList[j]-1) - 1;
						indexB = momId[indexList[j]-1] - 1;
					}
					auto partialsIdxA = partialsIdxMap[j].at(indexA);
					auto partialsIdxB = partialsIdxMap[j].at(indexB);
					for (int i = 0; i < psDim; i++) {
						partials[i][partialsIdxA] -= value * sum[i][j][0];
						partials[i][partialsIdxB] -= value * sum[i][j][1];
					}
				}
			});

	return;
}

void PSISuperCluster::computeDissociationPartialDerivatives(double* partials[5],
		const ReactionNetwork::PartialsIdxMap partialsIdxMap[5]) const {

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)

	// Visit all the dissociating pairs
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this,
			&partials, &partialsIdxMap](DissociationPairMap::value_type const& currMapItem) {
				auto& currPair = currMapItem.second;

				// Get the dissociating clusters
				auto const& cluster = currPair.first;
				// Compute the contribution from the dissociating cluster
				auto value = currPair.kConstant / (double) nTot;
				for (int j = 0; j < psDim; j++) {
					int index = 0;
					if (j == 0) {
						index = cluster.getId() - 1;
					}
					else {
						index = cluster.getMomentId(indexList[j]-1) - 1;
					}
					auto partialsIdx = partialsIdxMap[j].at(index);
					for (int i = 0; i < psDim; i++) {
						partials[i][partialsIdx] += value * currPair.coefs[j][i];
					}
				}
			});

	return;
}

void PSISuperCluster::computeEmissionPartialDerivatives(double* partials[5],
		const ReactionNetwork::PartialsIdxMap partialsIdxMap[5]) const {

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)

	// Visit all the emission pairs
	std::for_each(effEmissionList.begin(), effEmissionList.end(),
			[this,
			&partials, &partialsIdxMap](DissociationPairMap::value_type const& currMapItem) {
				auto& currPair = currMapItem.second;

				// Compute the contribution from the dissociating cluster
				auto value = currPair.kConstant / (double) nTot;
				for (int j = 0; j < psDim; j++) {
					int index = 0;
					if (j == 0) {
						index = id - 1;
					}
					else {
						index = momId[indexList[j]-1] - 1;
					}
					auto partialsIdx = partialsIdxMap[j].at(index);
					for (int i = 0; i < psDim; i++) {
						partials[i][partialsIdx] -= value * currPair.coefs[j][i];
					}
				}
			});

	return;
}

void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::ProductionCoefficientBase const& curr) const {

	os << "a[0-4][0-4][0-4]: ";
	for (const auto& curr2D : curr.coefs) {
		for (const auto& curr1D : curr2D) {
			std::copy(curr1D.begin(), curr1D.end(),
					std::ostream_iterator<double>(os, " "));
		}
	}
}

void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::SuperClusterDissociationPair const& currPair) const {

	os << "a[0-4][0-4]: ";
	for (const auto& curr1D : currPair.coefs) {
		std::copy(curr1D.begin(), curr1D.end(),
				std::ostream_iterator<double>(os, " "));
	}
}

void PSISuperCluster::outputCoefficientsTo(std::ostream& os) const {

	os << "id: " << id << '\n';
	os << "reacting: " << effReactingList.size() << '\n';
	std::for_each(effReactingList.begin(), effReactingList.end(),
			[this,&os](ProductionPairMap::value_type const& currMapItem) {
				auto const& currPair = currMapItem.second;
				os << "first: " << currPair.first.getId()
				<< "; second: " << currPair.second.getId()
				<< "; ";
				dumpCoefficients(os, currPair);
				os << '\n';
			});

	os << "combining: " << effCombiningList.size() << '\n';
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[this,&os](CombiningClusterMap::value_type const& currMapItem) {
				auto const& currComb = currMapItem.second;
				os << "other: " << currComb.first.getId()
				<< "; ";
				dumpCoefficients(os, currComb);
				os << '\n';
			});

	os << "dissociating: " << effDissociatingList.size() << '\n';
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this,&os](DissociationPairMap::value_type const& currMapItem) {
				auto const& currPair = currMapItem.second;
				os << "first: " << currPair.first.getId()
				<< "; second: " << currPair.second.getId()
				<< "; ";
				dumpCoefficients(os, currPair);
				os << '\n';
			});

	os << "emitting: " << effEmissionList.size() << '\n';
	std::for_each(effEmissionList.begin(), effEmissionList.end(),
			[this,&os](DissociationPairMap::value_type const& currMapItem) {
				auto const& currPair = currMapItem.second;
				os << "first: " << currPair.first.getId()
				<< "; second: " << currPair.second.getId()
				<< "; ";
				dumpCoefficients(os, currPair);
				os << '\n';
			});
}

