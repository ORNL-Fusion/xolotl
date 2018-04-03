// Includes
#include "PSISuperCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <MathUtils.h>
#include <xolotlPerf.h>

using namespace xolotlCore;

/**
 * The moment partials.
 */
std::vector<double> psiMomentPartials[4];

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
		for (int i = 0; i < 4; i++) {
			firstDistance[i + 1] = super.getDistance(b[i], i);
		}
	}
	if (prodPair.second.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(prodPair.second);
		for (int i = 0; i < 4; i++) {
			secondDistance[i + 1] = super.getDistance(b[i], i);
		}
	}
	double factor[5] = { };
	factor[0] = 1.0;
	for (int i = 0; i < 4; i++) {
		factor[i + 1] = getFactor(a[i], i);
	}
	// First is A, second is B, in A + B -> this
	for (int k = 0; k < 5; k++) {
		for (int j = 0; j < 5; j++) {
			for (int i = 0; i < 5; i++) {
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
					for (int i = 0; i < 4; i++) {
						firstDistance[i + 1] = super.getDistance(currPRI.b[i], i);
					}
				}
				if (prodPair.second.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster const&>(prodPair.second);
					for (int i = 0; i < 4; i++) {
						secondDistance[i + 1] = super.getDistance(currPRI.b[i], i);
					}
				}
				double factor[5] = {};
				factor[0] = 1.0;
				for (int i = 0; i < 4; i++) {
					factor[i + 1] = getFactor(currPRI.a[i], i);
				}
				// First is A, second is B, in A + B -> this
				for (int k = 0; k < 5; k++) {
					for (int j = 0; j < 5; j++) {
						for (int i = 0; i < 5; i++) {
							prodPair.coefs[i][j][k] += firstDistance[i] * secondDistance[j] * factor[k];
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
	for (int i = 0; i < 4; i++) {
		distance[i + 1] = getDistance(a[i], i);
		factor[i + 1] = getFactor(a[i], i);
	}

	// This is A, itBis is B, in A + B -> C
	for (int k = 0; k < 5; k++) {
		for (int j = 0; j < 5; j++) {
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
				for (int i = 0; i < 4; i++) {
					distance[i+1] = getDistance(currPRInfo.b[i], i);
					factor[i+1] = getFactor(currPRInfo.b[i], i);
				}

				// This is A, itBis is B, in A + B -> C
				for (int k = 0; k < 5; k++) {
					for (int j = 0; j < 5; j++) {
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
		for (int i = 0; i < 4; i++) {
			distance[i + 1] = super.getDistance(a[i], i);
		}
	}
	for (int i = 0; i < 4; i++) {
		factor[i + 1] = getFactor(b[i], i);
	}
	// A is the dissociating cluster
	for (int j = 0; j < 5; j++) {
		for (int i = 0; i < 5; i++) {
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
					for (int i = 0; i < 4; i++) {
						distance[i + 1] = super.getDistance(currPRI.a[i], i);
					}
				}
				for (int i = 0; i < 4; i++) {
					factor[i + 1] = getFactor(currPRI.b[i], i);
				}
				// A is the dissociating cluster
				for (int j = 0; j < 5; j++) {
					for (int i = 0; i < 5; i++) {
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
	for (int i = 0; i < 4; i++) {
		distance[i + 1] = getDistance(a[i], i);
		factor[i + 1] = getFactor(a[i], i);
	}
	// A is the dissociating cluster
	for (int j = 0; j < 5; j++) {
		for (int i = 0; i < 5; i++) {
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
				for (int i = 0; i < 4; i++) {
					distance[i+1] = getDistance(currPRI.a[i], i);
					factor[i+1] = getFactor(currPRI.a[i], i);
				}
				// A is the dissociating cluster
				for (int j = 0; j < 5; j++) {
					for (int i = 0; i < 5; i++) {
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
		reactionRadius += xolotlCore::tungstenLatticeConstant
				* pow((3.0 * (double) (std::get<3>(pair))) / xolotlCore::pi,
						(1.0 / 3.0)) * 0.5 / (double) nTot;

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

		// Add the concentration of each cluster in the group times its number of helium
		conc += getConcentration(heDistance, dDistance, tDistance, vDistance)
				* (double) std::get<3>(pair);
	}

	return conc;
}

double PSISuperCluster::getIntegratedVConcentration(int v) const {
	// Initial declarations
	int heIndex = 0, dIndex = 0, tIndex = 0;
	double heDistance = 0.0, dDistance = 0.0, tDistance = 0.0, vDistance = 0.0,
			conc = 0.0;

	// Loop on the widths
	for (int l = 0; l < sectionWidth[2]; l++) {
		// Compute the tritium index
		tIndex = (int) (numAtom[2] - (double) sectionWidth[2] / 2.0) + l + 1;
		for (int k = 0; k < sectionWidth[1]; k++) {
			// Compute the deuterium index
			dIndex = (int) (numAtom[1] - (double) sectionWidth[1] / 2.0) + k
					+ 1;
			for (int j = 0; j < sectionWidth[0]; j++) {
				// Compute the helium index
				heIndex = (int) (numAtom[0] - (double) sectionWidth[0] / 2.0)
						+ j + 1;

				// Check if this cluster exists
				if (heVList.find(std::make_tuple(heIndex, dIndex, tIndex, v))
						== heVList.end())
					continue;

				// Compute the distances
				heDistance = getDistance(heIndex, 0);
				dDistance = getDistance(dIndex, 1);
				tDistance = getDistance(tIndex, 2);
				vDistance = getDistance(v, 3);

				// Add the concentration of each cluster in the group times its number of helium
				conc += getConcentration(heDistance, dDistance, tDistance,
						vDistance);
			}
		}
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
				for (int i = 0; i < 4; i++) {
					setReactionConnectivity(currPair.first.getMomentId(i));
					setReactionConnectivity(currPair.second.getMomentId(i));
				}
			});

	// Visit all the combining pairs
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[this](CombiningClusterMap::value_type const& currMapItem) {
				// The cluster is connecting to the combining cluster
				auto const& currComb = currMapItem.second;
				setReactionConnectivity(currComb.first.getId());
				for (int i = 0; i < 4; i++) {
					setReactionConnectivity(currComb.first.getMomentId(i));
				}
			});

	// Loop over all the dissociating pairs
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this](DissociationPairMap::value_type const& currMapItem) {
				// The cluster is connecting to the combining cluster
				auto const& currPair = currMapItem.second;
				setDissociationConnectivity(currPair.first.getId());
				for (int i = 0; i < 4; i++) {
					setDissociationConnectivity(currPair.first.getMomentId(i));
				}
			});

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	// Initialize the partial vector for the moment
	int dof = network.getDOF();
	psiMomentPartials[0].resize(dof, 0.0);
	psiMomentPartials[1].resize(dof, 0.0);
	psiMomentPartials[2].resize(dof, 0.0);
	psiMomentPartials[3].resize(dof, 0.0);

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
				for (int i = 0; i < 4; i++) {
					lA[i+1] = dissociatingCluster.getMoment(i);
				}

				double sum[5] = {};
				for (int j = 0; j < 5; j++) {
					for (int i = 0; i < 5; i++) {
						sum[j] += currPair.coefs[i][j] * lA[i];
					}
				}
				// Update the flux
				auto value = currPair.kConstant / (double) nTot;
				flux += value * sum[0];
				// Compute the moment fluxes
				for (int i = 0; i < 4; i++) {
					momentFlux[i] += value * sum[i+1];
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
				for (int i = 0; i < 4; i++) {
					lA[i+1] = l1[i];
				}

				double sum[5] = {};
				for (int j = 0; j < 5; j++) {
					for (int i = 0; i < 5; i++) {
						sum[j] += currPair.coefs[i][j] * lA[i];
					}
				}
				// Update the flux
				auto value = currPair.kConstant / (double) nTot;
				flux += value * sum[0];
				// Compute the moment fluxes
				for (int i = 0; i < 4; i++) {
					momentFlux[i] -= value * sum[i+1];
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
				for (int i = 0; i < 4; i++) {
					lA[i+1] = firstReactant.getMoment(i);
					lB[i+1] = secondReactant.getMoment(i);
				}

				double sum[5] = {};
				for (int k = 0; k < 5; k++) {
					for (int j = 0; j < 5; j++) {
						for (int i = 0; i < 5; i++) {
							sum[k] += currPair.coefs[j][i][k] * lA[j] * lB[i];
						}
					}
				}

				// Update the flux
				auto value = currPair.kConstant / (double) nTot;
				flux += value * sum[0];
				// Compute the moment fluxes
				for (int i = 0; i < 4; i++) {
					momentFlux[i] += value * sum[i+1];
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
				for (int i = 0; i < 4; i++) {
					lA[i+1] = l1[i];
					lB[i+1] = combiningCluster.getMoment(i);
				}

				double sum[5] = {};
				for (int k = 0; k < 5; k++) {
					for (int j = 0; j < 5; j++) {
						for (int i = 0; i < 5; i++) {
							sum[k] += currComb.coefs[i][j][k] * lA[i] * lB[j];
						}
					}
				}
				// Update the flux
				auto value = currComb.kConstant / (double) nTot;
				flux += value * sum[0];
				// Compute the moment fluxes
				for (int i = 0; i < 4; i++) {
					momentFlux[i] -= value * sum[i+1];
				}
			});

	return flux;
}

void PSISuperCluster::getPartialDerivatives(
		std::vector<double> & partials) const {
	// Reinitialize the moment partial derivatives vector
	std::fill(psiMomentPartials[0].begin(), psiMomentPartials[0].end(), 0.0);
	std::fill(psiMomentPartials[1].begin(), psiMomentPartials[1].end(), 0.0);
	std::fill(psiMomentPartials[2].begin(), psiMomentPartials[2].end(), 0.0);
	std::fill(psiMomentPartials[3].begin(), psiMomentPartials[3].end(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return;
}

void PSISuperCluster::getProductionPartialDerivatives(
		std::vector<double> & partials) const {

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A

	// Loop over all the reacting pairs
	std::for_each(effReactingList.begin(), effReactingList.end(),
			[this,&partials](ProductionPairMap::value_type const& currMapItem) {

				auto const& currPair = currMapItem.second;
				// Get the two reacting clusters
				auto const& firstReactant = currPair.first;
				auto const& secondReactant = currPair.second;
				double lA[5] = {}, lB[5] = {};
				lA[0] = firstReactant.getConcentration();
				lB[0] = secondReactant.getConcentration();
				for (int i = 0; i < 4; i++) {
					lA[i+1] = firstReactant.getMoment(i);
					lB[i+1] = secondReactant.getMoment(i);
				}

				double sum[5][5][2] = {};
				for (int k = 0; k < 5; k++) {
					for (int j = 0; j < 5; j++) {
						for (int i = 0; i < 5; i++) {
							sum[k][j][0] += currPair.coefs[j][i][k] * lB[i];
							sum[k][j][1] += currPair.coefs[i][j][k] * lA[i];
						}
					}
				}

				// Compute the contribution from the first and second part of the reacting pair
				auto value = currPair.kConstant / (double) nTot;
				for (int j = 0; j < 5; j++) {
					int indexA = 0, indexB = 0;
					if (j == 0) {
						indexA = firstReactant.getId() - 1;
						indexB = secondReactant.getId() - 1;
					}
					else {
						indexA = firstReactant.getMomentId(j-1) - 1;
						indexB = secondReactant.getMomentId(j-1) - 1;
					}
					partials[indexA] += value * sum[0][j][0];
					partials[indexB] += value * sum[0][j][1];
					for (int i = 0; i < 4; i++) {
						psiMomentPartials[i][indexA] += value * sum[i+1][j][0];
						psiMomentPartials[i][indexB] += value * sum[i+1][j][1];
					}
				}
			});

	return;
}

void PSISuperCluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A

	// Visit all the combining clusters
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[this,&partials](CombiningClusterMap::value_type const& currMapItem) {
				// Get the combining clusters
				auto const& currComb = currMapItem.second;
				auto const& cluster = currComb.first;
				double lA[5] = {}, lB[5] = {};
				lA[0] = l0;
				lB[0] = cluster.getConcentration();
				for (int i = 0; i < 4; i++) {
					lA[i+1] = l1[i];
					lB[i+1] = cluster.getMoment(i);
				}

				double sum[5][5][2] = {};
				for (int k = 0; k < 5; k++) {
					for (int j = 0; j < 5; j++) {
						for (int i = 0; i < 5; i++) {
							sum[k][j][0] += currComb.coefs[i][j][k] * lA[i];
							sum[k][j][1] += currComb.coefs[j][i][k] * lB[i];
						}
					}
				}

				// Compute the contribution from the both clusters
				auto value = currComb.kConstant / (double) nTot;
				for (int j = 0; j < 5; j++) {
					int indexA = 0, indexB = 0;
					if (j == 0) {
						indexA = cluster.getId() - 1;
						indexB = id - 1;
					}
					else {
						indexA = cluster.getMomentId(j-1) - 1;
						indexB = momId[j-1] - 1;
					}
					partials[indexA] -= value * sum[0][j][0];
					partials[indexB] -= value * sum[0][j][1];
					for (int i = 0; i < 4; i++) {
						psiMomentPartials[i][indexA] -= value * sum[i+1][j][0];
						psiMomentPartials[i][indexB] -= value * sum[i+1][j][1];
					}
				}
			});

	return;
}

void PSISuperCluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)

	// Visit all the dissociating pairs
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this,&partials](DissociationPairMap::value_type const& currMapItem) {
				auto& currPair = currMapItem.second;

				// Get the dissociating clusters
				auto const& cluster = currPair.first;
				// Compute the contribution from the dissociating cluster
				auto value = currPair.kConstant / (double) nTot;
				for (int j = 0; j < 5; j++) {
					int index = 0;
					if (j == 0) {
						index = cluster.getId() - 1;
					}
					else {
						index = cluster.getMomentId(j-1) - 1;
					}
					partials[index] += value * currPair.coefs[j][0];
					for (int i = 0; i < 4; i++) {
						psiMomentPartials[i][index] += value * currPair.coefs[j][i];
					}
				}
			});

	return;
}

void PSISuperCluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)

	// Visit all the emission pairs
	std::for_each(effEmissionList.begin(), effEmissionList.end(),
			[this,&partials](DissociationPairMap::value_type const& currMapItem) {
				auto& currPair = currMapItem.second;

				// Compute the contribution from the dissociating cluster
				auto value = currPair.kConstant / (double) nTot;
				for (int j = 0; j < 5; j++) {
					int index = 0;
					if (j == 0) {
						index = id - 1;
					}
					else {
						index = momId[j-1] - 1;
					}
					partials[index] -= value * currPair.coefs[j][0];
					for (int i = 0; i < 4; i++) {
						psiMomentPartials[i][index] -= value * currPair.coefs[j][i];
					}
				}
			});

	return;
}

void PSISuperCluster::getMomentPartialDerivatives(
		std::vector<double> & partials, int axis) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = psiMomentPartials[axis][i];
	}

	return;
}

void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::ProductionCoefficientBase const& curr) const {

	os << "a[0-4][0-4][0-4]: ";
	for (int k = 0; k < 5; k++) {
		for (int j = 0; j < 5; j++) {
			for (int i = 0; i < 5; i++) {
				os << curr.coefs[k][j][i] << ' ';
			}
		}
	}
}

void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::SuperClusterDissociationPair const& currPair) const {

	os << "a[0-4][0-4]: ";

	for (int j = 0; j < 5; j++) {
		for (int i = 0; i < 5; i++) {
			os << currPair.coefs[j][i] << ' ';
		}
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

