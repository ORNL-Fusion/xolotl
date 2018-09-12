// Includes
#include <iterator>
#include "PSISuperCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <xolotlPerf.h>

using namespace xolotlCore;

PSISuperCluster::PSISuperCluster(double num[4], int _nTot, int width[4],
		int lower[4], int higher[4], IReactionNetwork& _network,
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

		// Set the boundaries
		bounds[i] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(lower[i]),
				static_cast<IReactant::SizeType>(higher[i] + 1));
	}

	// Update the composition map
	composition[toCompIdx(Species::He)] = (int) num[0];
	composition[toCompIdx(Species::D)] = (int) num[1];
	composition[toCompIdx(Species::T)] = (int) num[2];
	composition[toCompIdx(Species::V)] = (int) num[3];
	composition[toCompIdx(Species::I)] = _nTot;

	// Set the formation energy
	formationEnergy = 0.0; // It is set to 0.0 because we do not want the super clusters to undergo dissociation

	// Set the diffusion factor and the migration energy
	migrationEnergy = std::numeric_limits<double>::infinity();
	diffusionFactor = 0.0;

	// Set the typename appropriately
	type = ReactantType::PSISuper;

	// Check the shape of the cluster
	full = (sectionWidth[0] * sectionWidth[1] * sectionWidth[2]
			* sectionWidth[3] == nTot) ? true : false;

	return;
}

auto PSISuperCluster::addToEffReactingList(ProductionReaction& reaction) 
        -> ProductionPairList::iterator {

	// Check if we already know about the reaction.
    auto listit = effReactingList.end();
	auto rkey = std::make_pair(&(reaction.first), &(reaction.second));
	auto mapit = effReactingListMap.find(rkey);
    if (mapit != effReactingListMap.end()) {
        // We already knew about this reaction.
        listit = mapit->second;
    }
    else {
        // We did not already know about this reaction.
        // Add info about reaction to our list.
        effReactingList.emplace_back(reaction,
						static_cast<PSICluster&>(reaction.first),
						static_cast<PSICluster&>(reaction.second), psDim);
        listit = effReactingList.end();
        --listit;
    }
	assert(listit != effReactingList.end());
    return listit;
}

auto PSISuperCluster::addToEffCombiningList(ProductionReaction& reaction) 
        -> CombiningClusterList::iterator {

	// Determine which is the "other" cluster.
    auto& otherCluster = findOtherCluster(reaction);

    // Check if we already know about the reaction.
    auto listit = effCombiningList.end();
	auto rkey = &otherCluster;
	auto mapit = effCombiningListMap.find(rkey);
    if (mapit != effCombiningListMap.end()) {
        // We already knew about this reaction.
        listit = mapit->second;
    }
    else {
		// We did not already know about the reaction.
		// Add info about reaction to our list.
        effCombiningList.emplace_back(reaction,
						static_cast<PSICluster&>(otherCluster), psDim);
        listit = effCombiningList.end();
        --listit;
	}
	assert(listit != effCombiningList.end());
    return listit;
}

auto PSISuperCluster::addToEffDissociatingList(DissociationReaction& reaction)
        -> DissociationPairList::iterator {

	// Determine which is the other cluster.
    auto& emittedCluster = findOtherCluster(reaction);

	// Check if we already know about the reaction.
    auto listit = effDissociatingList.end();
	auto rkey = std::make_pair(&(reaction.dissociating), &emittedCluster);
	auto mapit = effDissociatingListMap.find(rkey);
    if (mapit != effDissociatingListMap.end()) {

        // We already knew about this reaction.
        listit = mapit->second;
    }
    else {

		// We did not already know about it.
		// Add info about reaction to our list.
        effDissociatingList.emplace_back(reaction,
						static_cast<PSICluster&>(reaction.dissociating),
						static_cast<PSICluster&>(emittedCluster), psDim);
        listit = effDissociatingList.end();
        --listit;
	}
	assert(listit != effDissociatingList.end());
    return listit;
}


auto PSISuperCluster::addToEffEmissionList(DissociationReaction& reaction)
        -> DissociationPairList::iterator {

	// Check if we already know about the reaction.
    auto listit = effEmissionList.end();
	auto rkey = std::make_pair(&(reaction.first), &(reaction.second));
	auto mapit = effEmissionListMap.find(rkey);
	if (mapit != effEmissionListMap.end()) {
        // We already knew about the reaction.
        listit = mapit->second;
    }
    else {

		// We did not already know about it.
		// Note that we emit from the two rectants according to the given
		// reaction.
        effEmissionList.emplace_back(reaction,
						static_cast<PSICluster&>(reaction.first),
						static_cast<PSICluster&>(reaction.second), psDim);
        listit = effEmissionList.end();
        --listit;
	}
	assert(listit != effEmissionList.end());
    return listit;
}




void PSISuperCluster::resultFrom(ProductionReaction& reaction, int a[4],
		int b[4]) {

    // Ensure we know about the reaction.
    auto listit = addToEffReactingList(reaction);
	auto& prodPair = *listit;

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
    auto listit = addToEffReactingList(reaction);
	auto& prodPair = *(listit);

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

void PSISuperCluster::resultFrom(ProductionReaction& reaction,
		IReactant& product) {

	// Check if we already know about the reaction.
    auto listit = addToEffReactingList(reaction);
	auto& prodPair = *(listit);

	auto const& superR1 = static_cast<PSICluster const&>(prodPair.first);
	auto const& superR2 = static_cast<PSICluster const&>(prodPair.second);
	auto const& superProd = static_cast<PSICluster const&>(product);

	// Check if an interstitial cluster is involved
	int iSize = 0;
	if (superR1.getType() == ReactantType::I) {
		iSize = superR1.getSize();
	} else if (superR2.getType() == ReactantType::I) {
		iSize = superR2.getSize();
	}

	// Loop on the different type of clusters in grouping
	int productLo[4] = { }, productHi[4] = { }, singleComp[4] = { }, r1Lo[4] =
			{ }, r1Hi[4] = { }, width[4] = { };
	int nOverlap = 1;
	for (int i = 1; i < 5; i++) {
		// Check the boundaries in all the directions
		auto const& bounds = superProd.getBounds(i - 1);
		productLo[i - 1] = *(bounds.begin()), productHi[i - 1] = *(bounds.end())
				- 1;

		if (prodPair.first.getType() == ReactantType::PSISuper) {
			auto const& r1Bounds = superR1.getBounds(i - 1);
			r1Lo[i - 1] = *(r1Bounds.begin()), r1Hi[i - 1] = *(r1Bounds.end())
					- 1;
			auto const& r2Bounds = superR2.getBounds(i - 1);
			singleComp[i - 1] = *(r2Bounds.begin());
		}

		if (prodPair.second.getType() == ReactantType::PSISuper) {
			auto const& r1Bounds = superR1.getBounds(i - 1);
			singleComp[i - 1] = *(r1Bounds.begin());
			auto const& r2Bounds = superR2.getBounds(i - 1);
			r1Lo[i - 1] = *(r2Bounds.begin()), r1Hi[i - 1] = *(r2Bounds.end())
					- 1;
		}

		// Special case for V and I
		if (i == 4)
			singleComp[i - 1] -= iSize;

		width[i - 1] = std::min(productHi[i - 1],
				r1Hi[i - 1] + singleComp[i - 1])
				- std::max(productLo[i - 1], r1Lo[i - 1] + singleComp[i - 1])
				+ 1;

		nOverlap *= width[i - 1];
	}

	// Compute the coefficients
	prodPair.coefs[0][0][0] += (double) nOverlap;
	for (int i = 1; i < psDim; i++) {
		prodPair.coefs[0][0][i] += ((double) nOverlap
				/ (dispersion[indexList[i] - 1] * (double) width[i - 1]))
				* firstOrderSum(
						std::max(productLo[i - 1],
								r1Lo[i - 1] + singleComp[i - 1]),
						std::min(productHi[i - 1],
								r1Hi[i - 1] + singleComp[i - 1]),
						numAtom[indexList[i] - 1]);

		double a = 0.0;
		if (r1Hi[i - 1] != r1Lo[i - 1]) {
			a = ((double) (nOverlap * 2)
					/ (double) ((r1Hi[i - 1] - r1Lo[i - 1]) * width[i - 1]))
					* firstOrderSum(
							std::max(productLo[i - 1] - singleComp[i - 1],
									r1Lo[i - 1]),
							std::min(productHi[i - 1] - singleComp[i - 1],
									r1Hi[i - 1]),
							(double) (r1Lo[i - 1] + r1Hi[i - 1]) / 2.0);

			prodPair.coefs[0][i][0] += prodPair.second.isMixed() * a;
			prodPair.coefs[i][0][0] += prodPair.first.isMixed() * a;

			a = ((double) (nOverlap * 2)
					/ ((double) ((r1Hi[i - 1] - r1Lo[i - 1]) * width[i - 1])
							* dispersion[indexList[i] - 1]))
					* secondOrderOffsetSum(
							std::max(productLo[i - 1],
									r1Lo[i - 1] + singleComp[i - 1]),
							std::min(productHi[i - 1],
									r1Hi[i - 1] + singleComp[i - 1]),
							numAtom[indexList[i] - 1],
							(double) (r1Lo[i - 1] + r1Hi[i - 1]) / 2.0,
							-singleComp[i - 1]);

			prodPair.coefs[0][i][i] += prodPair.second.isMixed() * a;
			prodPair.coefs[i][0][i] += prodPair.first.isMixed() * a;
		}

		for (int j = 1; j < psDim; j++) {
			if (i == j)
				continue;

			if (r1Hi[i - 1] != r1Lo[i - 1]) {
				a = ((double) (nOverlap * 2)
						/ ((double) ((r1Hi[i - 1] - r1Lo[i - 1]) * width[i - 1]
								* width[j - 1]) * dispersion[indexList[j] - 1]))
						* firstOrderSum(
								std::max(productLo[i - 1] - singleComp[i - 1],
										r1Lo[i - 1]),
								std::min(productHi[i - 1] - singleComp[i - 1],
										r1Hi[i - 1]),
								(double) (r1Lo[i - 1] + r1Hi[i - 1]) / 2.0)
						* firstOrderSum(
								std::max(productLo[j - 1],
										r1Lo[j - 1] + singleComp[j - 1]),
								std::min(productHi[j - 1],
										r1Hi[j - 1] + singleComp[j - 1]),
								numAtom[indexList[j] - 1]);

				prodPair.coefs[0][i][j] += prodPair.second.isMixed() * a;
				prodPair.coefs[i][0][j] += prodPair.first.isMixed() * a;
			}
		}
	}

	return;
}

void PSISuperCluster::resultFrom(ProductionReaction& reaction, double *coef) {

	// Check if we already know about the reaction.
    auto listit = addToEffReactingList(reaction);
	auto& prodPair = *(listit);

	// NB: newPair's reactants are same as reaction's.
	// So use newPair only from here on.
	// TODO Any way to enforce this beyond splitting it into two functions?

	// Update the coefficients
	int n = 0;
	for (int i = 0; i < psDim; i++) {
		for (int j = 0; j < psDim; j++) {
			for (int k = 0; k < psDim; k++) {
				prodPair.coefs[i][j][k] += coef[n];
				n++;
			}
		}
	}

	return;
}


void PSISuperCluster::participateIn(ProductionReaction& reaction, int a[4]) {

	// Ensure we know about the reaction.
    auto listit = addToEffCombiningList(reaction);
	auto& combCluster = *listit;

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

	// Ensure we know about the reaction.
    auto listit = addToEffCombiningList(reaction);
	auto& combCluster = *listit;

	// Update the coefficients
	std::for_each(pendingPRInfos.begin(), pendingPRInfos.end(),
			[this,&combCluster](const PendingProductionReactionInfo& currPRInfo) {
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

void PSISuperCluster::participateIn(ProductionReaction& reaction,
		IReactant& product) {

	// Ensure we know about the reaction.
    auto listit = addToEffCombiningList(reaction);
	auto& combCluster = *listit;

	auto const& superProd = static_cast<PSICluster const&>(product);

	// Check if an interstitial cluster is involved
    auto& otherCluster = findOtherCluster(reaction);
	int iSize = 0;
	if (otherCluster.getType() == ReactantType::I) {
		iSize = otherCluster.getSize();
	}

	// Loop on the different type of clusters in grouping
	int productLo[4] = { }, productHi[4] = { }, singleComp[4] = { }, r1Lo[4] =
			{ }, r1Hi[4] = { }, width[4] = { };
	int nOverlap = 1;
	for (int i = 1; i < 5; i++) {
		auto const& bounds = superProd.getBounds(i - 1);
		productLo[i - 1] = *(bounds.begin()), productHi[i - 1] = *(bounds.end())
				- 1;
		auto const& r1Bounds = getBounds(i - 1);
		r1Lo[i - 1] = *(r1Bounds.begin()), r1Hi[i - 1] = *(r1Bounds.end()) - 1;
		auto const& r2Bounds = otherCluster.getBounds(i - 1);
		singleComp[i - 1] = *(r2Bounds.begin());

		// Special case for V and I
		if (i == 4)
			singleComp[i - 1] -= iSize;

		width[i - 1] = std::min(productHi[i - 1],
				r1Hi[i - 1] + singleComp[i - 1])
				- std::max(productLo[i - 1], r1Lo[i - 1] + singleComp[i - 1])
				+ 1;

		nOverlap *= width[i - 1];
	}

	// Compute the coefficients
	combCluster.coefs[0][0][0] += nOverlap;
	for (int i = 1; i < psDim; i++) {
		combCluster.coefs[0][0][i] += ((double) nOverlap
				/ (dispersion[indexList[i] - 1] * (double) width[i - 1]))
				* firstOrderSum(
						std::max(productLo[i - 1] - singleComp[i - 1],
								r1Lo[i - 1]),
						std::min(productHi[i - 1] - singleComp[i - 1],
								r1Hi[i - 1]), numAtom[indexList[i] - 1]);

		if (sectionWidth[indexList[i] - 1] != 1)
			combCluster.coefs[i][0][0] += ((double) (nOverlap * 2)
					/ (double) ((sectionWidth[indexList[i] - 1] - 1)
							* width[i - 1]))
					* firstOrderSum(
							std::max(productLo[i - 1] - singleComp[i - 1],
									r1Lo[i - 1]),
							std::min(productHi[i - 1] - singleComp[i - 1],
									r1Hi[i - 1]), numAtom[indexList[i] - 1]);

		if (sectionWidth[indexList[i] - 1] != 1)
			combCluster.coefs[i][0][i] += ((double) (nOverlap * 2)
					/ ((double) ((sectionWidth[indexList[i] - 1] - 1)
							* width[i - 1]) * dispersion[indexList[i] - 1]))
					* secondOrderSum(
							std::max(productLo[i - 1] - singleComp[i - 1],
									r1Lo[i - 1]),
							std::min(productHi[i - 1] - singleComp[i - 1],
									r1Hi[i - 1]), numAtom[indexList[i] - 1]);

		for (int j = 1; j < psDim; j++) {
			if (i == j)
				continue;

			if (sectionWidth[indexList[i] - 1] != 1)
				combCluster.coefs[i][0][j] += ((double) (nOverlap * 2)
						/ ((double) ((sectionWidth[indexList[i] - 1] - 1)
								* width[i - 1] * width[j - 1])
								* dispersion[indexList[j] - 1]))
						* firstOrderSum(
								std::max(productLo[i - 1] - singleComp[i - 1],
										r1Lo[i - 1]),
								std::min(productHi[i - 1] - singleComp[i - 1],
										r1Hi[i - 1]), numAtom[indexList[i] - 1])
						* firstOrderSum(
								std::max(productLo[j - 1] - singleComp[j - 1],
										r1Lo[j - 1]),
								std::min(productHi[j - 1] - singleComp[j - 1],
										r1Hi[j - 1]),
								numAtom[indexList[j] - 1]);
		}
	}

	return;
}

void PSISuperCluster::participateIn(ProductionReaction& reaction,
		double *coef) {


	// Ensure we know about the reaction.
    auto listit = addToEffCombiningList(reaction);
	auto& combCluster = *listit;

	// Update the coefficients
	int n = 0;
	for (int i = 0; i < psDim; i++) {
		for (int j = 0; j < psDim; j++) {
			for (int k = 0; k < psDim; k++) {
				combCluster.coefs[i][j][k] += coef[n];
				n++;
			}
		}
	}

	return;
}

void PSISuperCluster::participateIn(DissociationReaction& reaction, int a[4],
		int b[4]) {

    // Ensure we know about the reaction.
    auto listit = addToEffDissociatingList(reaction);
	auto& dissPair = *listit;

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

    // Ensure we know about the reaction.
    auto listit = addToEffDissociatingList(reaction);
    auto& dissPair = *listit;

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

void PSISuperCluster::participateIn(DissociationReaction& reaction,
		IReactant& disso) {

    // Ensure we know about the reaction.
    auto listit = addToEffDissociatingList(reaction);
	auto& dissPair = *listit;

	auto const& superDisso = static_cast<PSICluster const&>(disso);

	// Check if an interstitial cluster is involved
    auto& emittedCluster = findOtherCluster(reaction);
	int iSize = 0;
	if (emittedCluster.getType() == ReactantType::I) {
		iSize = emittedCluster.getSize();
	}

	// Loop on the different type of clusters in grouping
	int dissoLo[4] = { }, dissoHi[4] = { }, singleComp[4] = { }, r1Lo[4] = { },
			r1Hi[4] = { }, width[4] = { };
	int nOverlap = 1;
	for (int i = 1; i < 5; i++) {
		auto const& bounds = superDisso.getBounds(i - 1);
		dissoLo[i - 1] = *(bounds.begin()), dissoHi[i - 1] = *(bounds.end()) - 1;
		auto const& r1Bounds = getBounds(i - 1);
		r1Lo[i - 1] = *(r1Bounds.begin()), r1Hi[i - 1] = *(r1Bounds.end()) - 1;
		auto const& r2Bounds = emittedCluster.getBounds(i - 1);
		singleComp[i - 1] = *(r2Bounds.begin());

		// Special case for V and I
		if (i == 4)
			singleComp[i - 1] -= iSize;

		width[i - 1] = std::min(dissoHi[i - 1], r1Hi[i - 1] + singleComp[i - 1])
				- std::max(dissoLo[i - 1], r1Lo[i - 1] + singleComp[i - 1]) + 1;

		nOverlap *= width[i - 1];
	}

	// Compute the coefficients
	dissPair.coefs[0][0] += nOverlap;
	for (int i = 1; i < psDim; i++) {
		dissPair.coefs[0][i] += ((double) nOverlap
				/ (dispersion[indexList[i] - 1] * (double) width[i - 1]))
				* firstOrderSum(
						std::max(dissoLo[i - 1] - singleComp[i - 1],
								r1Lo[i - 1]),
						std::min(dissoHi[i - 1] - singleComp[i - 1],
								r1Hi[i - 1]), numAtom[indexList[i] - 1]);

		if (dissoHi[i - 1] != dissoLo[i - 1]) {
			dissPair.coefs[i][0] +=
					((double) (2 * nOverlap)
							/ (double) ((dissoHi[i - 1] - dissoLo[i - 1])
									* width[i - 1]))
							* firstOrderSum(
									std::max(dissoLo[i - 1],
											r1Lo[i - 1] + singleComp[i - 1]),
									std::min(dissoHi[i - 1],
											r1Hi[i - 1] + singleComp[i - 1]),
									(double) (dissoLo[i - 1] + dissoHi[i - 1])
											/ 2.0);

			dissPair.coefs[i][i] += ((double) (2 * nOverlap)
					/ ((double) ((dissoHi[i - 1] - dissoLo[i - 1])
							* width[i - 1]) * dispersion[indexList[i] - 1]))
					* secondOrderOffsetSum(
							std::max(dissoLo[i - 1],
									r1Lo[i - 1] + singleComp[i - 1]),
							std::min(dissoHi[i - 1],
									r1Hi[i - 1] + singleComp[i - 1]),
							(double) (dissoLo[i - 1] + dissoHi[i - 1]) / 2.0,
							numAtom[indexList[i] - 1], -singleComp[i - 1]);
		}

		for (int j = 1; j < psDim; j++) {
			if (i == j)
				continue;

			if (dissoHi[i - 1] != dissoLo[i - 1])
				dissPair.coefs[i][j] += ((double) (nOverlap * 2)
						/ ((double) ((dissoHi[i - 1] - dissoLo[i - 1])
								* width[i - 1] * width[j - 1])
								* dispersion[indexList[j] - 1]))
						* firstOrderSum(
								std::max(dissoLo[i - 1],
										singleComp[i - 1] + r1Lo[i - 1]),
								std::min(dissoHi[i - 1],
										singleComp[i - 1] + r1Hi[i - 1]),
								(double) (dissoLo[i - 1] + dissoHi[i - 1])
										/ 2.0)
						* firstOrderSum(
								std::max(dissoLo[j - 1] - singleComp[j - 1],
										r1Lo[j - 1]),
								std::min(dissoHi[j - 1] - singleComp[j - 1],
										r1Hi[j - 1]),
								numAtom[indexList[j] - 1]);
		}
	}

	return;
}

void PSISuperCluster::participateIn(DissociationReaction& reaction,
		double *coef) {

    // Ensure we know about the reaction.
    auto listit = addToEffDissociatingList(reaction);
	auto& dissPair = *listit;

	// Update the coefficients
	int n = 0;
	for (int i = 0; i < psDim; i++) {
		for (int j = 0; j < psDim; j++) {
			dissPair.coefs[i][j] += coef[n];
			n++;
		}
	}

	return;
}


void PSISuperCluster::emitFrom(DissociationReaction& reaction, int a[4]) {

    // Ensure we know about the reaction.
    auto listit = addToEffEmissionList(reaction);
	auto& dissPair = *listit;

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

    // Ensure we know about the reaction.
    auto listit = addToEffEmissionList(reaction);
	auto& dissPair = *listit;


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

void PSISuperCluster::emitFrom(DissociationReaction& reaction,
		IReactant& disso) {

    // Ensure we know about the reaction.
    auto listit = addToEffEmissionList(reaction);
	auto& dissPair = *listit;


	auto const& superR1 = static_cast<PSICluster const&>(dissPair.first);
	auto const& superR2 = static_cast<PSICluster const&>(dissPair.second);
	auto const& superDisso = static_cast<PSICluster const&>(disso);

	// Check if an interstitial cluster is involved
	int iSize = 0;
	if (superR1.getType() == ReactantType::I) {
		iSize = superR1.getSize();
	} else if (superR2.getType() == ReactantType::I) {
		iSize = superR2.getSize();
	}

	// Loop on the different type of clusters in grouping
	int dissoLo[4] = { }, dissoHi[4] = { }, singleComp[4] = { }, r1Lo[4] = { },
			r1Hi[4] = { }, width[4] = { };
	int nOverlap = 1;
	for (int i = 1; i < 5; i++) {
		// Check the boundaries in all the directions
		auto const& bounds = superDisso.getBounds(i - 1);
		dissoLo[i - 1] = *(bounds.begin()), dissoHi[i - 1] = *(bounds.end()) - 1;

		if (dissPair.first.getType() == ReactantType::PSISuper) {
			auto const& r1Bounds = superR1.getBounds(i - 1);
			r1Lo[i - 1] = *(r1Bounds.begin()), r1Hi[i - 1] = *(r1Bounds.end())
					- 1;
			auto const& r2Bounds = superR2.getBounds(i - 1);
			singleComp[i - 1] = *(r2Bounds.begin());
		}

		if (dissPair.second.getType() == ReactantType::PSISuper) {
			auto const& r1Bounds = superR1.getBounds(i - 1);
			singleComp[i - 1] = *(r1Bounds.begin());
			auto const& r2Bounds = superR2.getBounds(i - 1);
			r1Lo[i - 1] = *(r2Bounds.begin()), r1Hi[i - 1] = *(r2Bounds.end())
					- 1;
		}

		// Special case for V and I
		if (i == 4)
			singleComp[i - 1] -= iSize;

		width[i - 1] = std::min(dissoHi[i - 1], r1Hi[i - 1] + singleComp[i - 1])
				- std::max(dissoLo[i - 1], r1Lo[i - 1] + singleComp[i - 1]) + 1;

		nOverlap *= width[i - 1];
	}

	// Compute the coefficients
	dissPair.coefs[0][0] += (double) nOverlap;
	for (int i = 1; i < psDim; i++) {
		dissPair.coefs[0][i] += ((double) nOverlap
				/ (dispersion[indexList[i] - 1] * (double) width[i - 1]))
				* firstOrderSum(
						std::max(dissoLo[i - 1],
								r1Lo[i - 1] + singleComp[i - 1]),
						std::min(dissoHi[i - 1],
								r1Hi[i - 1] + singleComp[i - 1]),
						numAtom[indexList[i] - 1]);

		if (sectionWidth[indexList[i] - 1] != 1) {
			dissPair.coefs[i][0] += ((double) (2 * nOverlap)
					/ (double) ((sectionWidth[indexList[i] - 1] - 1)
							* width[i - 1]))
					* firstOrderSum(
							std::max(dissoLo[i - 1],
									r1Lo[i - 1] + singleComp[i - 1]),
							std::min(dissoHi[i - 1],
									r1Hi[i - 1] + singleComp[i - 1]),
							numAtom[indexList[i] - 1]);

			dissPair.coefs[i][i] += ((double) (2 * nOverlap)
					/ ((double) ((sectionWidth[indexList[i] - 1] - 1)
							* width[i - 1]) * dispersion[indexList[i] - 1]))
					* secondOrderSum(
							std::max(dissoLo[i - 1],
									r1Lo[i - 1] + singleComp[i - 1]),
							std::min(dissoHi[i - 1],
									r1Hi[i - 1] + singleComp[i - 1]),
							numAtom[indexList[i] - 1]);
		}

		for (int j = 1; j < psDim; j++) {
			if (i == j)
				continue;

			if (sectionWidth[indexList[i] - 1] != 1)
				dissPair.coefs[i][j] += ((double) (2 * nOverlap)
						/ ((double) (width[i - 1] * width[j - 1]
								* (sectionWidth[indexList[i] - 1] - 1))
								* dispersion[indexList[j] - 1]))
						* firstOrderSum(
								std::max(dissoLo[i - 1],
										r1Lo[i - 1] + singleComp[i - 1]),
								std::min(dissoHi[i - 1],
										r1Hi[i - 1] + singleComp[i - 1]),
								numAtom[indexList[i] - 1])
						* firstOrderSum(
								std::max(dissoLo[j - 1],
										r1Lo[j - 1] + singleComp[j - 1]),
								std::min(dissoHi[j - 1],
										r1Hi[j - 1] + singleComp[j - 1]),
								numAtom[indexList[j] - 1]);
		}
	}

	return;
}

void PSISuperCluster::emitFrom(DissociationReaction& reaction, double *coef) {

    // Ensure we know about the reaction.
    auto listit = addToEffEmissionList(reaction);
	auto& dissPair = *listit;

	// Update the coefficients
	int n = 0;
	for (int i = 0; i < psDim; i++) {
		for (int j = 0; j < psDim; j++) {
			dissPair.coefs[i][j] += coef[n];
			n++;
		}
	}

	return;
}

void PSISuperCluster::setHeVVector(const HeVListType& vec) {

	// Copy the list of coordinates
	heVList = vec;

	// Initialize the dispersion sum
	double nSquare[4] = { };
	// Update the network map, compute the radius and dispersions
	for (auto const& pair : heVList) {
        constexpr auto tlcCubed = xolotlCore::tungstenLatticeConstant *
                                    xolotlCore::tungstenLatticeConstant *
                                    xolotlCore::tungstenLatticeConstant;
		double rad = (sqrt(3.0) / 4.0) * xolotlCore::tungstenLatticeConstant
				+ cbrt((3.0 * tlcCubed * std::get<3>(pair)) / (8.0 * xolotlCore::pi))
				- cbrt((3.0 * tlcCubed) / (8.0 * xolotlCore::pi));
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
			[this](ProductionPairList::value_type const& currPair) {
				// The cluster is connecting to both clusters in the pair
				setReactionConnectivity(currPair.first.getId());
				setReactionConnectivity(currPair.second.getId());
				for (int i = 1; i < psDim; i++) {
					setReactionConnectivity(currPair.first.getMomentId(indexList[i]-1));
					setReactionConnectivity(currPair.second.getMomentId(indexList[i]-1));
				}
			});

	// Visit all the combining pairs
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[this](CombiningClusterList::value_type const& currComb) {
				// The cluster is connecting to the combining cluster
				setReactionConnectivity(currComb.first.getId());
				for (int i = 1; i < psDim; i++) {
					setReactionConnectivity(currComb.first.getMomentId(indexList[i]-1));
				}
			});

	// Loop over all the dissociating pairs
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this](DissociationPairList::value_type const& currPair) {
				// The cluster is connecting to the combining cluster
				setDissociationConnectivity(currPair.first.getId());
				for (int i = 1; i < psDim; i++) {
					setDissociationConnectivity(currPair.first.getMomentId(indexList[i]-1));
				}
			});

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

    // We're done with the maps used to construct our effective reaction lists.
    // Release them to reclaim the memory.
    // TODO is this true going forward?  Any desire to add
    // more reactions?
    effReactingListMap.clear();
    effCombiningListMap.clear();
    effDissociatingListMap.clear();
    effEmissionListMap.clear();

	return;
}

double PSISuperCluster::getDissociationFlux(int xi) {
	// Initial declarations
	double flux = 0.0;

	// Sum over all the dissociating pairs
	// TODO consider using std::accumulate.  May also want to change side
	// effect of updating member variables heMomentFlux and
	// vMomentFlux here.
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this,&flux,&xi](DissociationPairList::value_type const& currPair) {

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
				auto value = currPair.reaction.kConstant[xi] / (double) nTot;
				flux += value * sum[0];
				// Compute the moment fluxes
				for (int i = 1; i < psDim; i++) {
					momentFlux[indexList[i]-1] += value * sum[i];
				}
			});

	// Return the flux
	return flux;
}

double PSISuperCluster::getEmissionFlux(int xi) {
	// Initial declarations
	double flux = 0.0;

	// Loop over all the emission pairs
	// TODO consider using std::accumulate.  May also want to change side
	// effect of updating member variables heMomentFlux and
	// vMomentFlux here.
	std::for_each(effEmissionList.begin(), effEmissionList.end(),
			[this,&flux,&xi](DissociationPairList::value_type const& currPair) {
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
				auto value = currPair.reaction.kConstant[xi] / (double) nTot;
				flux += value * sum[0];
				// Compute the moment fluxes
				for (int i = 1; i < psDim; i++) {
					momentFlux[indexList[i]-1] -= value * sum[i];
				}
			});

	return flux;
}

double PSISuperCluster::getProductionFlux(int xi) {
	// Local declarations
	double flux = 0.0;

	// Sum over all the reacting pairs
	// TODO consider using std::accumulate.  May also want to change side
	// effect of updating member variables heMomentFlux and
	// vMomentFlux here.
	std::for_each(effReactingList.begin(), effReactingList.end(),
			[this,&flux,&xi](ProductionPairList::value_type const& currPair) {

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
				auto value = currPair.reaction.kConstant[xi] / (double) nTot;
				flux += value * sum[0];
				// Compute the moment fluxes
				for (int i = 1; i < psDim; i++) {
					momentFlux[indexList[i]-1] += value * sum[i];
				}
			});

	// Return the production flux
	return flux;
}

double PSISuperCluster::getCombinationFlux(int xi) {
	// Local declarations
	double flux = 0.0;

	// Sum over all the combining clusters
	// TODO consider using std::accumulate.  May also want to change side
	// effect of updating member variables heMomentFlux and
	// vMomentFlux here.
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[this,&flux,&xi](CombiningClusterList::value_type const& currComb) {
				// Get the combining cluster
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
				auto value = currComb.reaction.kConstant[xi] / (double) nTot;
				flux += value * sum[0];
				// Compute the moment fluxes
				for (int i = 1; i < psDim; i++) {
					momentFlux[indexList[i]-1] -= value * sum[i];
				}
			});

	return flux;
}

void PSISuperCluster::computePartialDerivatives(double* partials[5],
		const ReactionNetwork::PartialsIdxMap partialsIdxMap[5], int i) const {

	// Get the partial derivatives for each reaction type
	computeProductionPartialDerivatives(partials, partialsIdxMap, i);
	computeCombinationPartialDerivatives(partials, partialsIdxMap, i);
	computeDissociationPartialDerivatives(partials, partialsIdxMap, i);
	computeEmissionPartialDerivatives(partials, partialsIdxMap, i);

	return;
}

void PSISuperCluster::computeProductionPartialDerivatives(double* partials[5],
		const ReactionNetwork::PartialsIdxMap partialsIdxMap[5], int xi) const {

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
			&partials, &partialsIdxMap,&xi](ProductionPairList::value_type const& currPair) {

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
				auto value = currPair.reaction.kConstant[xi] / (double) nTot;
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
		const ReactionNetwork::PartialsIdxMap partialsIdxMap[5], int xi) const {

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
			&partials, &partialsIdxMap,&xi](CombiningClusterList::value_type const& currComb) {
				// Get the combining clusters
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
				auto value = currComb.reaction.kConstant[xi] / (double) nTot;
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
		const ReactionNetwork::PartialsIdxMap partialsIdxMap[5], int xi) const {

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)

	// Visit all the dissociating pairs
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this,
			&partials, &partialsIdxMap,&xi](DissociationPairList::value_type const& currPair) {

				// Get the dissociating clusters
				auto const& cluster = currPair.first;
				// Compute the contribution from the dissociating cluster
				auto value = currPair.reaction.kConstant[xi] / (double) nTot;

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
		const ReactionNetwork::PartialsIdxMap partialsIdxMap[5], int xi) const {

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)

	// Visit all the emission pairs
	std::for_each(effEmissionList.begin(), effEmissionList.end(),
			[this,
			&partials, &partialsIdxMap,&xi](DissociationPairList::value_type const& currPair) {

				// Compute the contribution from the dissociating cluster
				auto value = currPair.reaction.kConstant[xi] / (double) nTot;
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

std::vector<std::vector<double> > PSISuperCluster::getProdVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the reacting pairs
	std::for_each(effReactingList.begin(), effReactingList.end(),
			[&toReturn,this](ProductionPairList::value_type const& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first.getId() - 1);
				tempVec.push_back(currPair.second.getId() - 1);
				for (int i = 0; i < psDim; i++) {
					for (int j = 0; j < psDim; j++) {
						for (int k = 0; k < psDim; k++) {
							tempVec.push_back(currPair.coefs[i][j][k]);
						}
					}
				}

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > PSISuperCluster::getCombVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the combining reactants
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[&toReturn,this](CombiningClusterList::value_type const& cc) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(cc.first.getId() - 1);
				for (int i = 0; i < psDim; i++) {
					for (int j = 0; j < psDim; j++) {
						for (int k = 0; k < psDim; k++) {
							tempVec.push_back(cc.coefs[i][j][k]);
						}
					}
				}

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > PSISuperCluster::getDissoVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the dissociating pairs
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[&toReturn,this](DissociationPairList::value_type const& currPair) {

				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first.getId() - 1);
				tempVec.push_back(currPair.second.getId() - 1);
				for (int i = 0; i < psDim; i++) {
					for (int j = 0; j < psDim; j++) {
						tempVec.push_back(currPair.coefs[i][j]);
					}
				}

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > PSISuperCluster::getEmitVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the emitting pairs
	std::for_each(effEmissionList.begin(), effEmissionList.end(),
			[&toReturn,this](DissociationPairList::value_type const& currPair) {

				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first.getId() - 1);
				tempVec.push_back(currPair.second.getId() - 1);
				for (int i = 0; i < psDim; i++) {
					for (int j = 0; j < psDim; j++) {
						tempVec.push_back(currPair.coefs[i][j]);
					}
				}

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::ProductionCoefficientBase const& curr) const {

	os << "a[0-4][0-4][0-4]: ";
	for (int k = 0; k < psDim; k++) {
		for (int j = 0; j < psDim; j++) {
			for (int i = 0; i < psDim; i++) {
				os << curr.coefs[k][j][i] << ' ';
			}
		}
	}
}

void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::SuperClusterDissociationPair const& currPair) const {

	os << "a[0-4][0-4]: ";
	for (int j = 0; j < psDim; j++) {
		for (int i = 0; i < psDim; i++) {
			os << currPair.coefs[j][i] << ' ';
		}
	}
}

void PSISuperCluster::outputCoefficientsTo(std::ostream& os) const {

	os << "name: " << name << '\n';
	os << "reacting: " << effReactingList.size() << '\n';
	std::for_each(effReactingList.begin(), effReactingList.end(),
			[this,&os](ProductionPairList::value_type const& currPair) {
				os << "first: " << currPair.first.getName()
				<< "; second: " << currPair.second.getName() << ";";
				dumpCoefficients(os, currPair);
				os << '\n';
			});

	os << "combining: " << effCombiningList.size() << '\n';
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[this,&os](CombiningClusterList::value_type const& currComb) {
				os << "other: " << currComb.first.getName() << ";";
				dumpCoefficients(os, currComb);
				os << '\n';
			});

	os << "dissociating: " << effDissociatingList.size() << '\n';
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this,&os](DissociationPairList::value_type const& currPair) {
				os << "first: " << currPair.first.getName()
				<< "; second: " << currPair.second.getName()
				<< "; ";
				dumpCoefficients(os, currPair);
				os << '\n';
			});

	os << "emitting: " << effEmissionList.size() << '\n';
	std::for_each(effEmissionList.begin(), effEmissionList.end(),
			[this,&os](DissociationPairList::value_type const& currPair) {
				os << "first: " << currPair.first.getName()
				<< "; second: " << currPair.second.getName()
				<< "; ";
				dumpCoefficients(os, currPair);
				os << '\n';
			});
}

