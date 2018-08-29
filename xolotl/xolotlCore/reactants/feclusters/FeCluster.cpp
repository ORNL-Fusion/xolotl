#include <algorithm>
#include <cassert>
#include "FeCluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

void FeCluster::resultFrom(ProductionReaction& reaction, int a[4], int b[4]) {

	// Add a cluster pair for the given reaction.
	reactingPairs.emplace_back(reaction,
			static_cast<FeCluster&>(reaction.first),
			static_cast<FeCluster&>(reaction.second));
	auto& newPair = reactingPairs.back();

	// NB: newPair's reactants are same as reaction's.
	// So use newPair only from here on.
	// TODO Any way to enforce this beyond splitting it into two functions?

	// Update the coefficients
	double firstHeDistance = 0.0, firstVDistance = 0.0, secondHeDistance = 0.0,
			secondVDistance = 0.0;
	if (newPair.first.getType() == ReactantType::FeSuper) {
		auto const& super = static_cast<FeCluster const&>(newPair.first);
		firstHeDistance = super.getHeDistance(b[0]);
		firstVDistance = super.getVDistance(b[1]);
	}
	if (newPair.second.getType() == ReactantType::FeSuper) {
		auto const& super = static_cast<FeCluster const&>(newPair.second);
		secondHeDistance = super.getHeDistance(b[0]);
		secondVDistance = super.getVDistance(b[1]);
	}
	newPair.a00 += 1.0;
	newPair.a10 += firstHeDistance;
	newPair.a20 += firstVDistance;
	newPair.a01 += secondHeDistance;
	newPair.a02 += secondVDistance;
	newPair.a11 += firstHeDistance * secondHeDistance;
	newPair.a12 += firstHeDistance * secondVDistance;
	newPair.a21 += firstVDistance * secondHeDistance;
	newPair.a22 += firstVDistance * secondVDistance;

	return;
}

void FeCluster::resultFrom(ProductionReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& prInfos) {

	// Add a cluster pair for the given reaction.
	reactingPairs.emplace_back(reaction,
			static_cast<FeCluster&>(reaction.first),
			static_cast<FeCluster&>(reaction.second));
	auto& newPair = reactingPairs.back();

	// NB: newPair's reactants are same as reaction's.
	// So use newPair only from here on.
	// TODO Any way to enforce this beyond splitting it into two functions?

	// Update the coefficients
	std::for_each(prInfos.begin(), prInfos.end(),
			[&newPair](const PendingProductionReactionInfo& currPRI) {

				// Use names that correspond to single version.
				int a = currPRI.a[0];
				int b = currPRI.a[1];
				int c = currPRI.b[0];
				int d = currPRI.b[1];

				double firstHeDistance = 0.0, firstVDistance = 0.0, secondHeDistance = 0.0,
				secondVDistance = 0.0;
				if (newPair.first.getType() == ReactantType::FeSuper) {
					auto const& super = static_cast<FeCluster const&>(newPair.first);
					firstHeDistance = super.getHeDistance(c);
					firstVDistance = super.getVDistance(d);
				}
				if (newPair.second.getType() == ReactantType::FeSuper) {
					auto const& super = static_cast<FeCluster const&>(newPair.second);
					secondHeDistance = super.getHeDistance(c);
					secondVDistance = super.getVDistance(d);
				}
				newPair.a00 += 1.0;
				newPair.a10 += firstHeDistance;
				newPair.a20 += firstVDistance;
				newPair.a01 += secondHeDistance;
				newPair.a02 += secondVDistance;
				newPair.a11 += firstHeDistance * secondHeDistance;
				newPair.a12 += firstHeDistance * secondVDistance;
				newPair.a21 += firstVDistance * secondHeDistance;
				newPair.a22 += firstVDistance * secondVDistance;
			});

	return;
}

void FeCluster::resultFrom(ProductionReaction& reaction, IReactant& product) {

	// Add a cluster pair for the given reaction.
	reactingPairs.emplace_back(reaction,
			static_cast<FeCluster&>(reaction.first),
			static_cast<FeCluster&>(reaction.second));
	auto& newPair = reactingPairs.back();

	// NB: newPair's reactants are same as reaction's.
	// So use newPair only from here on.
	// TODO Any way to enforce this beyond splitting it into two functions?

	// Values for grouping parameters
	int productLoHe = composition[toCompIdx(Species::He)], productHiHe =
			composition[toCompIdx(Species::He)], productLoV =
			composition[toCompIdx(Species::V)], productHiV =
			composition[toCompIdx(Species::V)], loHe = 0, hiHe = 0, loV = 0,
			hiV = 0, singleHeSize = 0, singleVSize = 0;

	if (newPair.first.getType() == ReactantType::FeSuper) {
		auto const& super = static_cast<FeCluster const&>(newPair.first);
		auto const& heBounds = super.getHeBounds();
		loHe = *(heBounds.begin());
		hiHe = *(heBounds.end()) - 1;
		auto const& vBounds = super.getVBounds();
		loV = *(vBounds.begin());
		hiV = *(vBounds.end()) - 1;
		auto singleComp = newPair.second.getComposition();
		singleHeSize = singleComp[toCompIdx(Species::He)];
		singleVSize = singleComp[toCompIdx(Species::V)]
				- singleComp[toCompIdx(Species::I)]; // can be < 0
	}
	if (newPair.second.getType() == ReactantType::FeSuper) {
		auto const& super = static_cast<FeCluster const&>(newPair.second);
		auto const& heBounds = super.getHeBounds();
		loHe = *(heBounds.begin());
		hiHe = *(heBounds.end()) - 1;
		auto const& vBounds = super.getVBounds();
		loV = *(vBounds.begin());
		hiV = *(vBounds.end()) - 1;
		auto singleComp = newPair.first.getComposition();
		singleHeSize = singleComp[toCompIdx(Species::He)];
		singleVSize = singleComp[toCompIdx(Species::V)]
				- singleComp[toCompIdx(Species::I)]; // can be < 0
	}

	int heWidth = std::min(productHiHe, hiHe + singleHeSize)
			- std::max(productLoHe, loHe + singleHeSize) + 1;
	int vWidth = std::min(productHiV, hiV + singleVSize)
			- std::max(productLoV, loV + singleVSize) + 1;

	newPair.a00 = heWidth * vWidth;

	newPair.a01 = ((double) (2 * vWidth) / (double) (hiHe - loHe))
			* firstOrderSum(std::max(productLoHe - singleHeSize, loHe),
					std::min(productHiHe - singleHeSize, hiHe),
					(double) (loHe + hiHe) / 2.0);

	newPair.a02 = ((double) (2 * heWidth) / (double) (hiV - loV))
			* firstOrderSum(std::max(productLoV - singleVSize, loV),
					std::min(productHiV - singleVSize, hiV),
					(double) (loV + hiV) / 2.0);

	return;
}

void FeCluster::resultFrom(ProductionReaction& reaction, double *coef) {

	// Add a cluster pair for the given reaction.
	reactingPairs.emplace_back(reaction,
			static_cast<FeCluster&>(reaction.first),
			static_cast<FeCluster&>(reaction.second));
	auto& newPair = reactingPairs.back();

	// Update the coefficients
	newPair.a00 += coef[0];
	newPair.a10 += coef[1];
	newPair.a20 += coef[2];
	newPair.a01 += coef[3];
	newPair.a11 += coef[4];
	newPair.a21 += coef[5];
	newPair.a02 += coef[6];
	newPair.a12 += coef[7];
	newPair.a22 += coef[8];

	return;
}

void FeCluster::participateIn(ProductionReaction& reaction, int a[4]) {
	// Look for the other cluster
	auto& otherCluster = static_cast<FeCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if the reaction was already added
	std::vector<CombiningCluster>::reverse_iterator it;
	for (it = combiningReactants.rbegin(); it != combiningReactants.rend();
			++it) {
		if (&otherCluster == &(it->combining)) {
			break;
		}
	}
	if (it == combiningReactants.rend()) {
		// We did not already know about this combination.
		// Note that we combine with the other cluster in this reaction.
		combiningReactants.emplace_back(reaction, otherCluster);
		it = combiningReactants.rbegin();
	}

	// Update the coefficients
	double heDistance = 0.0, vDistance = 0.0;
	if (otherCluster.getType() == ReactantType::FeSuper) {
		heDistance = otherCluster.getHeDistance(a[0]);
		vDistance = otherCluster.getVDistance(a[1]);
	}
	(*it).a0 += 1.0;
	(*it).a1 += heDistance;
	(*it).a2 += vDistance;

	return;
}

void FeCluster::participateIn(ProductionReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& prInfos) {
	// Look for the other cluster
	auto& otherCluster = static_cast<FeCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if the reaction was already added
	std::vector<CombiningCluster>::reverse_iterator it;
	for (it = combiningReactants.rbegin(); it != combiningReactants.rend();
			++it) {
		if (&otherCluster == &(it->combining)) {
			break;
		}
	}
	if (it == combiningReactants.rend()) {
		// We did not already know about this combination.
		// Note that we combine with the other cluster in this reaction.
		combiningReactants.emplace_back(reaction, otherCluster);
		it = combiningReactants.rbegin();
	}
	assert(it != combiningReactants.rend());
	auto& combCluster = *it;

	// Update the coefficients
	std::for_each(prInfos.begin(), prInfos.end(),
			[this,&otherCluster,&combCluster](const PendingProductionReactionInfo& currPRInfo) {

				// Use names corresponding to the single-item version.
				int a = currPRInfo.b[0];
				int b = currPRInfo.b[1];

				double heDistance = 0.0, vDistance = 0.0;
				if (otherCluster.getType() == ReactantType::FeSuper) {
					heDistance = otherCluster.getHeDistance(a);
					vDistance = otherCluster.getVDistance(b);
				}
				combCluster.a0 += 1.0;
				combCluster.a1 += heDistance;
				combCluster.a2 += vDistance;
			});

	return;
}

void FeCluster::participateIn(ProductionReaction& reaction, IReactant& prod) {

	// Look for the other cluster
	auto& otherCluster = static_cast<FeCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if the reaction was already added
	std::vector<CombiningCluster>::reverse_iterator it;
	for (it = combiningReactants.rbegin(); it != combiningReactants.rend();
			++it) {
		if (&otherCluster == &(it->combining)) {
			break;
		}
	}
	if (it == combiningReactants.rend()) {
		// We did not already know about this combination.
		// Note that we combine with the other cluster in this reaction.
		combiningReactants.emplace_back(reaction, otherCluster);
		it = combiningReactants.rbegin();
	}
	assert(it != combiningReactants.rend());
	auto& combCluster = *it;

	// Values for grouping parameters
	int singleHeSize = composition[toCompIdx(Species::He)];
	int singleVSize = composition[toCompIdx(Species::V)]
			- composition[toCompIdx(Species::I)];		// can be < 0

	auto const& other = static_cast<FeCluster const&>(otherCluster);
	auto const& heBounds = other.getHeBounds();
	int loHe = *(heBounds.begin());
	int hiHe = *(heBounds.end()) - 1;
	auto const& vBounds = other.getVBounds();
	int loV = *(vBounds.begin());
	int hiV = *(vBounds.end()) - 1;

	auto const& product = static_cast<FeCluster const&>(prod);
	auto const& productHeBounds = product.getHeBounds();
	int productLoHe = *(productHeBounds.begin());
	int productHiHe = *(productHeBounds.end()) - 1;
	auto const& productVBounds = product.getVBounds();
	int productLoV = *(productVBounds.begin());
	int productHiV = *(productVBounds.end()) - 1;

	int heWidth = std::min(productHiHe, hiHe + singleHeSize)
			- std::max(productLoHe, loHe + singleHeSize) + 1;
	int vWidth = std::min(productHiV, hiV + singleVSize)
			- std::max(productLoV, loV + singleVSize) + 1;

	combCluster.a0 += heWidth * vWidth;

	combCluster.a1 += ((double) (2 * vWidth) / (double) (hiHe - loHe))
			* firstOrderSum(std::max(productLoHe - singleHeSize, loHe),
					std::min(productHiHe - singleHeSize, hiHe),
					(double) (loHe + hiHe) / 2.0);

	combCluster.a2 += ((double) (2 * heWidth) / (double) (hiV - loV))
			* firstOrderSum(std::max(productLoV - singleVSize, loV),
					std::min(productHiV - singleVSize, hiV),
					(double) (loV + hiV) / 2.0);

	return;
}

void FeCluster::participateIn(ProductionReaction& reaction, double *coef) {
	// Look for the other cluster
	auto& otherCluster = static_cast<FeCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if the reaction was already added
	std::vector<CombiningCluster>::reverse_iterator it;
	for (it = combiningReactants.rbegin(); it != combiningReactants.rend();
			++it) {
		if (&otherCluster == &(it->combining)) {
			break;
		}
	}
	if (it == combiningReactants.rend()) {
		// We did not already know about this combination.
		// Note that we combine with the other cluster in this reaction.
		combiningReactants.emplace_back(reaction, otherCluster);
		it = combiningReactants.rbegin();
	}

	// Update the coefficients
	(*it).a0 += coef[0];
	(*it).a1 += coef[1];
	(*it).a2 += coef[2];

	return;
}

void FeCluster::participateIn(DissociationReaction& reaction, int a[4],
		int b[4]) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<FeCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if the reaction was already added
	auto it =
			std::find_if(dissociatingPairs.rbegin(), dissociatingPairs.rend(),
					[&reaction,&emittedCluster](const ClusterPair& currPair) {
						return (&(reaction.dissociating) == &static_cast<FeCluster&>(currPair.first)) and
						(&emittedCluster == &static_cast<FeCluster&>(currPair.second));

					});
	if (it == dissociatingPairs.rend()) {

		// We did not already know about it.

		// Add the pair of them where it is important that the
		// dissociating cluster is the first one
		dissociatingPairs.emplace_back(reaction,
				static_cast<FeCluster&>(reaction.dissociating),
				static_cast<FeCluster&>(emittedCluster));
		it = dissociatingPairs.rbegin();
	}

	// Update the coefficients
	double firstHeDistance = 0.0, firstVDistance = 0.0;
	if (reaction.dissociating.getType() == ReactantType::FeSuper) {
		auto const& super = static_cast<FeCluster&>(reaction.dissociating);
		firstHeDistance = super.getHeDistance(a[0]);
		firstVDistance = super.getVDistance(a[1]);
	}
	(*it).a00 += 1.0;
	(*it).a10 += firstHeDistance;
	(*it).a20 += firstVDistance;

	return;
}

void FeCluster::participateIn(DissociationReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& prInfos) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<FeCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if the reaction was already added
	auto it =
			std::find_if(dissociatingPairs.rbegin(), dissociatingPairs.rend(),
					[&reaction,&emittedCluster](const ClusterPair& currPair) {
						return (&(reaction.dissociating) == &static_cast<FeCluster&>(currPair.first)) and
						(&emittedCluster == &static_cast<FeCluster&>(currPair.second));

					});
	if (it == dissociatingPairs.rend()) {

		// We did not already know about it.

		// Add the pair of them where it is important that the
		// dissociating cluster is the first one
		dissociatingPairs.emplace_back(reaction,
				static_cast<FeCluster&>(reaction.dissociating),
				static_cast<FeCluster&>(emittedCluster));
		it = dissociatingPairs.rbegin();
	}
	assert(it != dissociatingPairs.rend());
	auto& currPair = *it;

	// Update the coefficients
	std::for_each(prInfos.begin(), prInfos.end(),
			[&currPair,&reaction](const PendingProductionReactionInfo& currPRI) {

				// Use names corresponding to the single-item version.
				int a = currPRI.a[0];
				int b = currPRI.a[1];

				double firstHeDistance = 0.0, firstVDistance = 0.0;
				if (reaction.dissociating.getType() == ReactantType::FeSuper) {
					auto const& super = static_cast<FeCluster&>(reaction.dissociating);
					firstHeDistance = super.getHeDistance(a);
					firstVDistance = super.getVDistance(b);
				}
				currPair.a00 += 1.0;
				currPair.a10 += firstHeDistance;
				currPair.a20 += firstVDistance;
			});

	return;
}

void FeCluster::participateIn(DissociationReaction& reaction,
		IReactant& disso) {

	// Look for the other cluster
	auto& emittedCluster = static_cast<FeCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if the reaction was already added
	auto it =
			std::find_if(dissociatingPairs.rbegin(), dissociatingPairs.rend(),
					[&reaction,&emittedCluster](const ClusterPair& currPair) {
						return (&(reaction.dissociating) == &static_cast<FeCluster&>(currPair.first)) and
						(&emittedCluster == &static_cast<FeCluster&>(currPair.second));

					});
	if (it == dissociatingPairs.rend()) {

		// We did not already know about it.

		// Add the pair of them where it is important that the
		// dissociating cluster is the first one
		dissociatingPairs.emplace_back(reaction,
				static_cast<FeCluster&>(reaction.dissociating),
				static_cast<FeCluster&>(emittedCluster));
		it = dissociatingPairs.rbegin();
	}
	assert(it != dissociatingPairs.rend());
	auto& currPair = *it;

	// Values for grouping parameters
	int singleHeSize = composition[toCompIdx(Species::He)];
	int singleVSize = composition[toCompIdx(Species::V)]
			- composition[toCompIdx(Species::I)]; // can be < 0

	auto const& other = static_cast<FeCluster const&>(emittedCluster);
	auto const& heBounds = other.getHeBounds();
	int loHe = *(heBounds.begin());
	int hiHe = *(heBounds.end()) - 1;
	auto const& vBounds = other.getVBounds();
	int loV = *(vBounds.begin());
	int hiV = *(vBounds.end()) - 1;

	auto const& dissoCluster = static_cast<FeCluster const&>(disso);
	auto const& dissoHeBounds = dissoCluster.getHeBounds();
	int dissoLoHe = *(dissoHeBounds.begin());
	int dissoHiHe = *(dissoHeBounds.end()) - 1;
	auto const& dissoVBounds = dissoCluster.getVBounds();
	int dissoLoV = *(dissoVBounds.begin());
	int dissoHiV = *(dissoVBounds.end()) - 1;

	int heWidth = std::min(dissoHiHe, hiHe + singleHeSize)
			- std::max(dissoLoHe, loHe + singleHeSize) + 1;
	int vWidth = std::min(dissoHiV, hiV + singleVSize)
			- std::max(dissoLoV, loV + singleVSize) + 1;

	currPair.a00 = heWidth * vWidth;

	currPair.a10 = ((double) (2 * vWidth) / (double) (dissoHiHe - dissoLoHe))
			* firstOrderSum(std::max(dissoLoHe, loHe + singleHeSize),
					std::min(dissoHiHe, hiHe + singleHeSize),
					(dissoLoHe + dissoHiHe) / 2.0);

	currPair.a20 = ((double) (2 * heWidth) / (double) (dissoHiV - dissoLoV))
			* firstOrderSum(std::max(dissoLoV, loV + singleVSize),
					std::min(dissoHiV, hiV + singleVSize),
					(dissoLoV + dissoHiV) / 2.0);

	if (disso.getType() != ReactantType::FeSuper) {
		currPair.a10 = 0.0, currPair.a20 = 0.0;
	}

	return;
}

void FeCluster::participateIn(DissociationReaction& reaction, double *coef) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<FeCluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if the reaction was already added
	auto it =
			std::find_if(dissociatingPairs.rbegin(), dissociatingPairs.rend(),
					[&reaction,&emittedCluster](const ClusterPair& currPair) {
						return (&(reaction.dissociating) == &static_cast<FeCluster&>(currPair.first)) and
						(&emittedCluster == &static_cast<FeCluster&>(currPair.second));

					});
	if (it == dissociatingPairs.rend()) {

		// We did not already know about it.

		// Add the pair of them where it is important that the
		// dissociating cluster is the first one
		dissociatingPairs.emplace_back(reaction,
				static_cast<FeCluster&>(reaction.dissociating),
				static_cast<FeCluster&>(emittedCluster));
		it = dissociatingPairs.rbegin();
	}

	// Update the coefficients
	(*it).a00 += coef[0];
	(*it).a10 += coef[1];
	(*it).a20 += coef[2];

	return;
}

void FeCluster::emitFrom(DissociationReaction& reaction, int a[4]) {

	// Note that we emit from the reaction's reactants according to
	// the given reaction.
	// TODO do we need to check to see whether we already know about
	// this reaction?
	emissionPairs.emplace_back(reaction,
			static_cast<FeCluster&>(reaction.first),
			static_cast<FeCluster&>(reaction.second));
	auto& dissPair = emissionPairs.back();

	// Count the number of reactions
	dissPair.a00 += 1.0;

	return;
}

void FeCluster::emitFrom(DissociationReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& prInfos) {

	// Note that we emit from the reaction's reactants according to
	// the given reaction.
	// TODO do we need to check to see whether we already know about
	// this reaction?
	emissionPairs.emplace_back(reaction,
			static_cast<FeCluster&>(reaction.first),
			static_cast<FeCluster&>(reaction.second));
	auto& dissPair = emissionPairs.back();

	// Update the coefficients
	std::for_each(prInfos.begin(), prInfos.end(),
			[&dissPair](const PendingProductionReactionInfo& currPRI) {
				// Update the coefficients
				dissPair.a00 += 1.0;
			});

	return;
}

void FeCluster::emitFrom(DissociationReaction& reaction, IReactant& disso) {

	// Note that we emit from the reaction's reactants according to
	// the given reaction.
	// TODO do we need to check to see whether we already know about
	// this reaction?
	emissionPairs.emplace_back(reaction,
			static_cast<FeCluster&>(reaction.first),
			static_cast<FeCluster&>(reaction.second));
	auto& dissPair = emissionPairs.back();

	// Values for grouping parameters
	int productLoHe = composition[toCompIdx(Species::He)], productHiHe =
			composition[toCompIdx(Species::He)], productLoV =
			composition[toCompIdx(Species::V)], productHiV =
			composition[toCompIdx(Species::V)], loHe = 0, hiHe = 0, loV = 0,
			hiV = 0, singleHeSize = 0, singleVSize = 0;

	if (dissPair.first.getType() == ReactantType::FeSuper) {
		auto const& super = static_cast<FeCluster const&>(dissPair.first);
		auto const& heBounds = super.getHeBounds();
		loHe = *(heBounds.begin());
		hiHe = *(heBounds.end()) - 1;
		auto const& vBounds = super.getVBounds();
		loV = *(vBounds.begin());
		hiV = *(vBounds.end()) - 1;
		auto singleComp = dissPair.second.getComposition();
		singleHeSize = singleComp[toCompIdx(Species::He)];
		singleVSize = singleComp[toCompIdx(Species::V)]
				- singleComp[toCompIdx(Species::I)]; // can be < 0
	}
	if (dissPair.second.getType() == ReactantType::FeSuper) {
		auto const& super = static_cast<FeCluster const&>(dissPair.second);
		auto const& heBounds = super.getHeBounds();
		loHe = *(heBounds.begin());
		hiHe = *(heBounds.end()) - 1;
		auto const& vBounds = super.getVBounds();
		loV = *(vBounds.begin());
		hiV = *(vBounds.end()) - 1;
		auto singleComp = dissPair.first.getComposition();
		singleHeSize = singleComp[toCompIdx(Species::He)];
		singleVSize = singleComp[toCompIdx(Species::V)]
				- singleComp[toCompIdx(Species::I)]; // can be < 0
	}

	int heWidth = std::min(productHiHe, hiHe + singleHeSize)
			- std::max(productLoHe, loHe + singleHeSize) + 1;
	int vWidth = std::min(productHiV, hiV + singleVSize)
			- std::max(productLoV, loV + singleVSize) + 1;

	dissPair.a00 = heWidth * vWidth;

	return;
}

void FeCluster::emitFrom(DissociationReaction& reaction, double *coef) {

	// Note that we emit from the reaction's reactants according to
	// the given reaction.
	// TODO do we need to check to see whether we already know about
	// this reaction?
	emissionPairs.emplace_back(reaction,
			static_cast<FeCluster&>(reaction.first),
			static_cast<FeCluster&>(reaction.second));
	auto& dissPair = emissionPairs.back();

	// Count the number of reactions
	dissPair.a00 += coef[0];

	return;
}

static std::vector<int> getFullConnectivityVector(std::set<int> connectivitySet,
		int size) {
	// Create a vector of zeroes with size equal to the network size
	std::vector<int> connectivity(size);

	// Set the value of the connectivity array to one for each element that is
	// in the set.
	for (auto it = connectivitySet.begin(); it != connectivitySet.end(); ++it) {
		connectivity[*it - 1] = 1;
	}

	return connectivity;
}

std::vector<int> FeCluster::getReactionConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(reactionConnectivitySet, network.getDOF());
}

std::vector<int> FeCluster::getDissociationConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(dissociationConnectivitySet,
			network.getDOF());
}

void FeCluster::resetConnectivities() {
	// Shrink the arrays to save some space. (About 10% or so.)
	reactingPairs.shrink_to_fit();
	combiningReactants.shrink_to_fit();
	dissociatingPairs.shrink_to_fit();
	emissionPairs.shrink_to_fit();

	// Clear both sets
	reactionConnectivitySet.clear();
	dissociationConnectivitySet.clear();

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);
	setDissociationConnectivity(id);
	setReactionConnectivity(momId[0]);
	setDissociationConnectivity(momId[0]);
	setReactionConnectivity(momId[1]);
	setDissociationConnectivity(momId[1]);

	// Loop on the effective reacting pairs
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[this](const ClusterPair& currPair) {
				// The cluster is connecting to both clusters in the pair
				setReactionConnectivity(currPair.first.id);
				setReactionConnectivity(currPair.second.id);
				setReactionConnectivity(currPair.first.momId[0]);
				setReactionConnectivity(currPair.second.momId[0]);
				setReactionConnectivity(currPair.first.momId[1]);
				setReactionConnectivity(currPair.second.momId[1]);
			});

	// Loop on the effective combining reactants
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[this](const CombiningCluster& cc) {
				// The cluster is connecting to the combining cluster
				setReactionConnectivity(cc.combining.id);
				setReactionConnectivity(cc.combining.momId[0]);
				setReactionConnectivity(cc.combining.momId[1]);
			});

	// Loop on the effective dissociating pairs
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[this](const ClusterPair& currPair) {
				// The cluster is connecting to the dissociating cluster which
				// is the first one by definition
				setDissociationConnectivity(currPair.first.id);
				setDissociationConnectivity(currPair.first.momId[0]);
				setDissociationConnectivity(currPair.first.momId[1]);
			});

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	return;
}

void FeCluster::updateFromNetwork() {

	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	return;
}

double FeCluster::getDissociationFlux(int xi) const {

	// Sum dissociation flux over all our dissociating clusters.
	double flux = std::accumulate(dissociatingPairs.begin(),
			dissociatingPairs.end(), 0.0,
			[&xi](double running, const ClusterPair& currPair) {
				auto const& dissCluster = currPair.first;
				double l0A = dissCluster.getConcentration(0.0, 0.0);
				double lHeA = dissCluster.getHeMoment();
				double lVA = dissCluster.getVMoment();

				// Calculate the Dissociation flux
				return running +
				(currPair.reaction.kConstant[xi] *
						(currPair.a00 * l0A +
								currPair.a10 * lHeA +
								currPair.a20 * lVA));
			});

	// Return the flux
	return flux;
}

double FeCluster::getEmissionFlux(int xi) const {

	// Sum rate constants from all emission pair reactions.
	double flux = std::accumulate(emissionPairs.begin(), emissionPairs.end(),
			0.0, [&xi](double running, const ClusterPair& currPair) {
				return running + currPair.reaction.kConstant[xi] * currPair.a00;
			});

	return flux * concentration;
}

double FeCluster::getProductionFlux(int xi) const {

	// Sum production flux over all reacting pairs.
	double flux = std::accumulate(reactingPairs.begin(), reactingPairs.end(),
			0.0, [&xi](double running, const ClusterPair& currPair) {

				// Get the two reacting clusters
			auto const& firstReactant = currPair.first;
			auto const& secondReactant = currPair.second;
			double l0A = firstReactant.getConcentration(0.0, 0.0);
			double l0B = secondReactant.getConcentration(0.0, 0.0);
			double lHeA = firstReactant.getHeMoment();
			double lHeB = secondReactant.getHeMoment();
			double lVA = firstReactant.getVMoment();
			double lVB = secondReactant.getVMoment();
			// Update the flux
			return running + currPair.reaction.kConstant[xi] *
			(currPair.a00 * l0A * l0B + currPair.a01 * l0A * lHeB +
					currPair.a02 * l0A * lVB + currPair.a10 * lHeA * l0B +
					currPair.a11 * lHeA * lHeB + currPair.a12 * lHeA * lVB +
					currPair.a20 * lVA * l0B + currPair.a21 * lVA * lHeB +
					currPair.a22 * lVA * lVB);
		});

	// Return the production flux
	return flux;
}

double FeCluster::getCombinationFlux(int xi) const {

	// Sum combination flux over all clusters that combine with us.
	double flux = std::accumulate(combiningReactants.begin(),
			combiningReactants.end(), 0.0,
			[&xi](double running, const CombiningCluster& cc) {

				// Get the cluster that combines with this one
				auto const& combiningCluster = cc.combining;
				double l0B = combiningCluster.getConcentration(0.0, 0.0);
				double lHeB = combiningCluster.getHeMoment();
				double lVB = combiningCluster.getVMoment();
				// Calculate the combination flux
				return running + (cc.reaction.kConstant[xi] *
						(cc.a0 * l0B + cc.a1 * lHeB + cc.a2 * lVB));

			});

	return flux * concentration;
}

std::vector<double> FeCluster::getPartialDerivatives(int i) const {
	// Local Declarations
	std::vector<double> partials(network.getDOF(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials, i);
	getCombinationPartialDerivatives(partials, i);
	getDissociationPartialDerivatives(partials, i);
	getEmissionPartialDerivatives(partials, i);

	return partials;
}

void FeCluster::getPartialDerivatives(std::vector<double> & partials,
		int i) const {
	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials, i);
	getCombinationPartialDerivatives(partials, i);
	getDissociationPartialDerivatives(partials, i);
	getEmissionPartialDerivatives(partials, i);

	return;
}

void FeCluster::getProductionPartialDerivatives(std::vector<double> & partials,
		int xi) const {

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[&partials,&xi](const ClusterPair& currPair) {
				// Get the two reacting clusters
				auto const& firstReactant = currPair.first;
				auto const& secondReactant = currPair.second;
				double l0A = firstReactant.getConcentration(0.0, 0.0);
				double l0B = secondReactant.getConcentration(0.0, 0.0);
				double lHeA = firstReactant.getHeMoment();
				double lHeB = secondReactant.getHeMoment();
				double lVA = firstReactant.getVMoment();
				double lVB = secondReactant.getVMoment();

				// Compute contribution from the first part of the reacting pair
				double value = currPair.reaction.kConstant[xi];

				partials[firstReactant.id - 1] += value *
				(currPair.a00 * l0B + currPair.a01 * lHeB + currPair.a02 * lVB);

				partials[firstReactant.momId[0] - 1] += value *
				(currPair.a10 * l0B + currPair.a11 * lHeB + currPair.a12 * lVB);

				partials[firstReactant.momId[1] - 1] += value *
				(currPair.a20 * l0B + currPair.a21 * lHeB + currPair.a22 * lVB);

				// Compute contribution from the second part of the reacting pair
				partials[secondReactant.id - 1] += value *
				(currPair.a00 * l0A + currPair.a10 * lHeA + currPair.a20 * lVA);

				partials[secondReactant.momId[0] - 1] += value *
				(currPair.a01 * l0A + currPair.a11 * lHeA + currPair.a21 * lVA);

				partials[secondReactant.momId[1] - 1] += value *
				(currPair.a02 * l0A + currPair.a12 * lHeA + currPair.a22 * lVA);
			});

	return;
}

void FeCluster::getCombinationPartialDerivatives(std::vector<double> & partials,
		int xi) const {

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[this,&partials,&xi](const CombiningCluster& cc) {
				auto const& cluster = cc.combining;
				double l0B = cluster.getConcentration(0.0, 0.0);
				double lHeB = cluster.getHeMoment();
				double lVB = cluster.getVMoment();

				// Remember that the flux due to combinations is OUTGOING (-=)!
				// Compute the contribution from this cluster
				partials[id - 1] -= cc.reaction.kConstant[xi]
				* (cc.a0 * l0B + cc.a1 * lHeB + cc.a2 * lVB);
				// Compute the contribution from the combining cluster
				double value = cc.reaction.kConstant[xi] * concentration;
				partials[cluster.id - 1] -= value * cc.a0;
				partials[cluster.momId[0] - 1] -= value * cc.a1;
				partials[cluster.momId[1] - 1] -= value * cc.a2;
			});

	return;
}

void FeCluster::getDissociationPartialDerivatives(
		std::vector<double> & partials, int xi) const {

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[&partials,&xi](const ClusterPair& currPair) {
				// Get the dissociating cluster
				auto const& cluster = currPair.first;
				double value = currPair.reaction.kConstant[xi];
				partials[cluster.id - 1] += value * currPair.a00;
				partials[cluster.momId[0] - 1] += value * currPair.a10;
				partials[cluster.momId[1] - 1] += value * currPair.a20;
			});

	return;
}

void FeCluster::getEmissionPartialDerivatives(std::vector<double> & partials,
		int xi) const {

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)
	double outgoingFlux = std::accumulate(emissionPairs.begin(),
			emissionPairs.end(), 0.0,
			[xi](double running, const ClusterPair& currPair) {
				return running + currPair.reaction.kConstant[xi] * currPair.a00;
			});
	partials[id - 1] -= outgoingFlux;

	return;
}

double FeCluster::getLeftSideRate(int i) const {

	// Sum rate constant-concentration product over combining reactants.
	double combiningRateTotal = std::accumulate(combiningReactants.begin(),
			combiningReactants.end(), 0.0,
			[&i](double running, const CombiningCluster& cc) {
				return running +
				(cc.reaction.kConstant[i] * cc.combining.concentration);
			});

	// Sum rate constants over all emission pair reactions.
	double emissionRateTotal = std::accumulate(emissionPairs.begin(),
			emissionPairs.end(), 0.0,
			[&i](double running, const ClusterPair& currPair) {
				return running + currPair.reaction.kConstant[i] * currPair.a00;
			});

	return combiningRateTotal + emissionRateTotal;
}

std::vector<std::vector<double> > FeCluster::getProdVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the reacting pairs
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[&toReturn](const ClusterPair& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first.getId() - 1);
				tempVec.push_back(currPair.second.getId() - 1);
				tempVec.push_back(currPair.a00);
				tempVec.push_back(currPair.a10);
				tempVec.push_back(currPair.a20);
				tempVec.push_back(currPair.a01);
				tempVec.push_back(currPair.a11);
				tempVec.push_back(currPair.a21);
				tempVec.push_back(currPair.a02);
				tempVec.push_back(currPair.a12);
				tempVec.push_back(currPair.a22);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > FeCluster::getCombVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the combining reactants
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[&toReturn](const CombiningCluster& cc) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(cc.combining.getId() - 1);
				tempVec.push_back(cc.a0);
				tempVec.push_back(cc.a1);
				tempVec.push_back(cc.a2);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > FeCluster::getDissoVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the dissociating pairs
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[&toReturn](const ClusterPair& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first.getId() - 1);
				tempVec.push_back(currPair.second.getId() - 1);
				tempVec.push_back(currPair.a00);
				tempVec.push_back(currPair.a10);
				tempVec.push_back(currPair.a20);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > FeCluster::getEmitVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the emitting pairs
	std::for_each(emissionPairs.begin(), emissionPairs.end(),
			[&toReturn](const ClusterPair& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first.getId() - 1);
				tempVec.push_back(currPair.second.getId() - 1);
				tempVec.push_back(currPair.a00);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<int> FeCluster::getConnectivity() const {
	int connectivityLength = network.getDOF();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);
	auto reactionConnVector = getReactionConnectivity();
	auto dissociationConnVector = getDissociationConnectivity();

	// The reaction and dissociation vectors must have a length equal to the
	// number of clusters
	if (reactionConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The reaction vector is an incorrect length");
	}
	if (dissociationConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The dissociation vector is an incorrect length");
	}

	// Merge the two vectors such that the final vector contains
	// a 1 at a position if either of the connectivity arrays
	// have a 1
	for (int i = 0; i < connectivityLength; i++) {
		// Consider each connectivity array only if its type is enabled
		connectivity[i] = reactionConnVector[i] || dissociationConnVector[i];
	}

	return connectivity;
}

void FeCluster::dumpCoefficients(std::ostream& os,
		FeCluster::ClusterPair const& curr) const {

	os << "a[0-2][0-2]: " << ' ' << curr.a00 << ' ' << curr.a01 << ' '
			<< curr.a02 << ' ' << curr.a10 << ' ' << curr.a11 << ' ' << curr.a12
			<< ' ' << curr.a20 << ' ' << curr.a21 << ' ' << curr.a22;
}

void FeCluster::dumpCoefficients(std::ostream& os,
		FeCluster::CombiningCluster const& curr) const {

	os << "a[0-2]: " << curr.a0 << ' ' << curr.a1 << ' ' << curr.a2;
}

void FeCluster::outputCoefficientsTo(std::ostream& os) const {

	os << "id: " << id << '\n';
	os << "reacting: " << reactingPairs.size() << '\n';
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[this,&os](ClusterPair const& currPair) {
				os << "first: " << currPair.first.getId()
				<< "; second: " << currPair.second.getId()
				<< "; ";
				dumpCoefficients(os, currPair);
				os << '\n';
			});

	os << "combining: " << combiningReactants.size() << '\n';
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[this,&os](CombiningCluster const& currCluster) {
				os << "other: " << currCluster.combining.getId()
				<< "; ";
				dumpCoefficients(os, currCluster);
				os << '\n';
			});

	os << "dissociating: " << dissociatingPairs.size() << '\n';
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[this,&os](ClusterPair const& currPair) {
				os << "first: " << currPair.first.getId()
				<< "; second: " << currPair.second.getId()
				<< "; ";
				dumpCoefficients(os, currPair);
				os << '\n';
			});

	os << "emitting: " << emissionPairs.size() << '\n';
	std::for_each(emissionPairs.begin(), emissionPairs.end(),
			[this,&os](ClusterPair const& currPair) {
				os << "first: " << currPair.first.getId()
				<< "; second: " << currPair.second.getId()
				<< "; ";
				dumpCoefficients(os, currPair);
				os << '\n';
			});
}

