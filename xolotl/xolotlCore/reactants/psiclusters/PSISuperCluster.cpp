// Includes
#include "PSISuperCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <MathUtils.h>
#include <xolotlPerf.h>

using namespace xolotlCore;

/**
 * The helium momentum partials.
 */
std::vector<double> psiHeMomentPartials;

/**
 * The vacancy momentum partials.
 */
std::vector<double> psiVMomentPartials;

PSISuperCluster::PSISuperCluster(double _numHe, double _numV, int _nTot,
		int heWidth, int vWidth, IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(_network, registry, buildName(_numHe, _numV)), numHe(_numHe), numV(
				_numV), nTot(_nTot), heBounds(0, 0), vBounds(0, 0), l0(0.0), l1He(
				0.0), l1V(0.0), dispersionHe(0.0), dispersionV(0.0), heMomentumFlux(
				0.0), vMomentumFlux(0.0) {
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = (int) (numHe + numV);

	// Update the composition map
	composition[toCompIdx(Species::He)] = (int) numHe;
	composition[toCompIdx(Species::V)] = (int) numV;

	// Set the width
	sectionHeWidth = heWidth;
	sectionVWidth = vWidth;

	// Set the formation energy
	formationEnergy = 0.0; // It is set to 0.0 because we do not want the super clusters to undergo dissociation

	// Set the diffusion factor and the migration energy
	migrationEnergy = std::numeric_limits<double>::infinity();
	diffusionFactor = 0.0;

	// Set the typename appropriately
	type = ReactantType::PSISuper;

	return;
}

void PSISuperCluster::resultFrom(ProductionReaction& reaction, int a, int b,
		int c, int d) {

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
	double firstHeDistance = 0.0, firstVDistance = 0.0, secondHeDistance = 0.0,
			secondVDistance = 0.0;
	if (prodPair.first.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(prodPair.first);
		firstHeDistance = super.getHeDistance(c);
		firstVDistance = super.getVDistance(d);
	}
	if (prodPair.second.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(prodPair.second);
		secondHeDistance = super.getHeDistance(c);
		secondVDistance = super.getVDistance(d);
	}
	double heFactor = (double) (a - numHe) / dispersionHe;
	double vFactor = (double) (b - numV) / dispersionV;
	// First is A, second is B, in A + B -> this
	prodPair.a000 += 1.0;
	prodPair.a001 += heFactor;
	prodPair.a002 += vFactor;
	prodPair.a100 += firstHeDistance;
	prodPair.a101 += firstHeDistance * heFactor;
	prodPair.a102 += firstHeDistance * vFactor;
	prodPair.a200 += firstVDistance;
	prodPair.a201 += firstVDistance * heFactor;
	prodPair.a202 += firstVDistance * vFactor;
	prodPair.a010 += secondHeDistance;
	prodPair.a011 += secondHeDistance * heFactor;
	prodPair.a012 += secondHeDistance * vFactor;
	prodPair.a020 += secondVDistance;
	prodPair.a021 += secondVDistance * heFactor;
	prodPair.a022 += secondVDistance * vFactor;
	prodPair.a110 += firstHeDistance * secondHeDistance;
	prodPair.a111 += firstHeDistance * secondHeDistance * heFactor;
	prodPair.a112 += firstHeDistance * secondHeDistance * vFactor;
	prodPair.a120 += firstHeDistance * secondVDistance;
	prodPair.a121 += firstHeDistance * secondVDistance * heFactor;
	prodPair.a122 += firstHeDistance * secondVDistance * vFactor;
	prodPair.a210 += firstVDistance * secondHeDistance;
	prodPair.a211 += firstVDistance * secondHeDistance * heFactor;
	prodPair.a212 += firstVDistance * secondHeDistance * vFactor;
	prodPair.a220 += firstVDistance * secondVDistance;
	prodPair.a221 += firstVDistance * secondVDistance * heFactor;
	prodPair.a222 += firstVDistance * secondVDistance * vFactor;

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

				// Use names corresponding to those in single version.
				int a = currPRI.numHe;
				int b = currPRI.numV;
				int c = currPRI.i;
				int d = currPRI.j;

				double firstHeDistance = 0.0, firstVDistance = 0.0, secondHeDistance = 0.0,
				secondVDistance = 0.0;
				if (prodPair.first.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster const&>(prodPair.first);
					firstHeDistance = super.getHeDistance(c);
					firstVDistance = super.getVDistance(d);
				}
				if (prodPair.second.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster const&>(prodPair.second);
					secondHeDistance = super.getHeDistance(c);
					secondVDistance = super.getVDistance(d);
				}
				double heFactor = (double) (a - numHe) / dispersionHe;
				double vFactor = (double) (b - numV) / dispersionV;
				// First is A, second is B, in A + B -> this
				prodPair.a000 += 1.0;
				prodPair.a001 += heFactor;
				prodPair.a002 += vFactor;
				prodPair.a100 += firstHeDistance;
				prodPair.a101 += firstHeDistance * heFactor;
				prodPair.a102 += firstHeDistance * vFactor;
				prodPair.a200 += firstVDistance;
				prodPair.a201 += firstVDistance * heFactor;
				prodPair.a202 += firstVDistance * vFactor;
				prodPair.a010 += secondHeDistance;
				prodPair.a011 += secondHeDistance * heFactor;
				prodPair.a012 += secondHeDistance * vFactor;
				prodPair.a020 += secondVDistance;
				prodPair.a021 += secondVDistance * heFactor;
				prodPair.a022 += secondVDistance * vFactor;
				prodPair.a110 += firstHeDistance * secondHeDistance;
				prodPair.a111 += firstHeDistance * secondHeDistance * heFactor;
				prodPair.a112 += firstHeDistance * secondHeDistance * vFactor;
				prodPair.a120 += firstHeDistance * secondVDistance;
				prodPair.a121 += firstHeDistance * secondVDistance * heFactor;
				prodPair.a122 += firstHeDistance * secondVDistance * vFactor;
				prodPair.a210 += firstVDistance * secondHeDistance;
				prodPair.a211 += firstVDistance * secondHeDistance * heFactor;
				prodPair.a212 += firstVDistance * secondHeDistance * vFactor;
				prodPair.a220 += firstVDistance * secondVDistance;
				prodPair.a221 += firstVDistance * secondVDistance * heFactor;
				prodPair.a222 += firstVDistance * secondVDistance * vFactor;
			});

	return;
}

void PSISuperCluster::participateIn(ProductionReaction& reaction, int a,
		int b) {

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
	double heDistance = getHeDistance(a);
	double heFactor = (double) (a - numHe) / dispersionHe;
	double vDistance = getVDistance(b);
	double vFactor = (double) (b - numV) / dispersionV;
	// This is A, itBis is B, in A + B -> C
	combCluster.a000 += 1.0;
	combCluster.a001 += heFactor;
	combCluster.a002 += vFactor;
	combCluster.a100 += heDistance;
	combCluster.a101 += heDistance * heFactor;
	combCluster.a102 += heDistance * vFactor;
	combCluster.a200 += vDistance;
	combCluster.a201 += vDistance * heFactor;
	combCluster.a202 += vDistance * vFactor;

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
			[this,&combCluster](const PendingProductionReactionInfo& currPRInfo) {

				// Use names corresponding to the single-item version.
				int a = currPRInfo.i;
				int b = currPRInfo.j;

				double heDistance = getHeDistance(a);
				double heFactor = (double) (a - numHe) / dispersionHe;
				double vDistance = getVDistance(b);
				double vFactor = (double) (b - numV) / dispersionV;
				// This is A, itBis is B, in A + B -> C
				combCluster.a000 += 1.0;
				combCluster.a001 += heFactor;
				combCluster.a002 += vFactor;
				combCluster.a100 += heDistance;
				combCluster.a101 += heDistance * heFactor;
				combCluster.a102 += heDistance * vFactor;
				combCluster.a200 += vDistance;
				combCluster.a201 += vDistance * heFactor;
				combCluster.a202 += vDistance * vFactor;
			});

	return;
}

void PSISuperCluster::participateIn(DissociationReaction& reaction, int a,
		int b, int c, int d) {

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
	double firstHeDistance = 0.0, firstVDistance = 0.0;
	if (reaction.dissociating.getType() == ReactantType::PSISuper) {
		auto const& super =
				static_cast<PSICluster const&>(reaction.dissociating);
		firstHeDistance = super.getHeDistance(a);
		firstVDistance = super.getVDistance(b);
	}
	double heFactor = (double) (c - numHe) / dispersionHe;
	double vFactor = (double) (d - numV) / dispersionV;

	// A is the dissociating cluster
	dissPair.a00 += 1.0;
	dissPair.a01 += heFactor;
	dissPair.a02 += vFactor;
	dissPair.a10 += firstHeDistance;
	dissPair.a11 += firstHeDistance * heFactor;
	dissPair.a12 += firstHeDistance * vFactor;
	dissPair.a20 += firstVDistance;
	dissPair.a21 += firstVDistance * heFactor;
	dissPair.a22 += firstVDistance * vFactor;

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

				// Use names corresponding to the single-item version.
				int a = currPRI.numHe;
				int b = currPRI.numV;
				int c = currPRI.i;
				int d = currPRI.j;

				double firstHeDistance = 0.0, firstVDistance = 0.0;
				if (reaction.dissociating.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster const&>(reaction.dissociating);
					firstHeDistance = super.getHeDistance(a);
					firstVDistance = super.getVDistance(b);
				}
				double heFactor = (double) (c - numHe) / dispersionHe;
				double vFactor = (double) (d - numV) / dispersionV;

				// A is the dissociating cluster
				dissPair.a00 += 1.0;
				dissPair.a01 += heFactor;
				dissPair.a02 += vFactor;
				dissPair.a10 += firstHeDistance;
				dissPair.a11 += firstHeDistance * heFactor;
				dissPair.a12 += firstHeDistance * vFactor;
				dissPair.a20 += firstVDistance;
				dissPair.a21 += firstVDistance * heFactor;
				dissPair.a22 += firstVDistance * vFactor;
			});

	return;
}

void PSISuperCluster::emitFrom(DissociationReaction& reaction, int a, int b,
		int c, int d) {

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

	// Update the coeeficients
	double heDistance = getHeDistance(a);
	double heFactor = (double) (a - numHe) / dispersionHe;
	double vDistance = getVDistance(b);
	double vFactor = (double) (b - numV) / dispersionV;
	// A is the dissociating cluster
	dissPair.a00 += 1.0;
	dissPair.a01 += heFactor;
	dissPair.a02 += vFactor;
	dissPair.a10 += heDistance;
	dissPair.a11 += heDistance * heFactor;
	dissPair.a12 += heDistance * vFactor;
	dissPair.a20 += vDistance;
	dissPair.a21 += vDistance * heFactor;
	dissPair.a22 += vDistance * vFactor;

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

	// Update the coeeficients
	std::for_each(prInfos.begin(), prInfos.end(),
			[this,&dissPair](const PendingProductionReactionInfo& currPRI) {

				// Use same names as used in single version.
				int a = currPRI.numHe;
				int b = currPRI.numV;

				double heDistance = getHeDistance(a);
				double heFactor = (double) (a - numHe) / dispersionHe;
				double vDistance = getVDistance(b);
				double vFactor = (double) (b - numV) / dispersionV;
				// A is the dissociating cluster
				dissPair.a00 += 1.0;
				dissPair.a01 += heFactor;
				dissPair.a02 += vFactor;
				dissPair.a10 += heDistance;
				dissPair.a11 += heDistance * heFactor;
				dissPair.a12 += heDistance * vFactor;
				dissPair.a20 += vDistance;
				dissPair.a21 += vDistance * heFactor;
				dissPair.a22 += vDistance * vFactor;
			});

	return;
}

void PSISuperCluster::setHeVVector(std::set<std::pair<int, int> > vec) {
	// Copy the list of coordinates
	heVList = vec;

	// Initialize the dispersion sum
	double nHeSquare = 0.0, nVSquare = 0.0;
	// Update the network map, compute the radius and dispersions
	for (auto it = heVList.begin(); it != heVList.end(); it++) {
		reactionRadius += xolotlCore::tungstenLatticeConstant
				* pow((3.0 * (double) ((*it).second)) / xolotlCore::pi,
						(1.0 / 3.0)) * 0.5 / (double) nTot;

		// Compute nSquare for the dispersion
		nHeSquare += (double) (*it).first * (*it).first;
		nVSquare += (double) (*it).second * (*it).second;
	}

	// Compute the dispersions
	if (sectionHeWidth == 1)
		dispersionHe = 1.0;
	else
		dispersionHe = 2.0 * (nHeSquare - (numHe * (double) nTot * numHe))
				/ ((double) (nTot * (sectionHeWidth - 1)));

	if (sectionVWidth == 1)
		dispersionV = 1.0;
	else
		dispersionV = 2.0 * (nVSquare - (numV * (double) nTot * numV))
				/ ((double) (nTot * (sectionVWidth - 1)));

	// Set the boundaries
	heBounds = IntegerRange<IReactant::SizeType>(
			static_cast<IReactant::SizeType>((numHe
					- (double) sectionHeWidth / 2.0) + 1),
			static_cast<IReactant::SizeType>((numHe
					- (double) sectionHeWidth / 2.0) + sectionHeWidth) + 1);
	vBounds = IntegerRange<IReactant::SizeType>(
			static_cast<IReactant::SizeType>((numV
					- (double) sectionVWidth / 2.0) + 1),
			static_cast<IReactant::SizeType>((numV
					- (double) sectionVWidth / 2.0) + sectionVWidth) + 1);

	return;
}

double PSISuperCluster::getTotalConcentration() const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
	for (auto const& pair : heVList) {
		// Compute the distances
		heDistance = getHeDistance(pair.first);
		vDistance = getVDistance(pair.second);

		// Add the concentration of each cluster in the group times its number of helium
		conc += getConcentration(heDistance, vDistance);
	}

	return conc;
}

double PSISuperCluster::getTotalHeliumConcentration() const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
	for (auto const& pair : heVList) {
		// Compute the distances
		heDistance = getHeDistance(pair.first);
		vDistance = getVDistance(pair.second);

		// Add the concentration of each cluster in the group times its number of helium
		conc += getConcentration(heDistance, vDistance) * (double) pair.first;
	}

	return conc;
}

double PSISuperCluster::getTotalVacancyConcentration() const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
	for (auto const& pair : heVList) {
		// Compute the distances
		heDistance = getHeDistance(pair.first);
		vDistance = getVDistance(pair.second);

		// Add the concentration of each cluster in the group times its number of vacancies
		conc += getConcentration(heDistance, vDistance) * (double) pair.second;
	}

	return conc;
}

double PSISuperCluster::getIntegratedVConcentration(int v) const {
	// Initial declarations
	int heIndex = 0;
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the helium width
	for (int j = 0; j < sectionHeWidth; j++) {
		// Compute the helium index
		heIndex = (int) (numHe - (double) sectionHeWidth / 2.0) + j + 1;

		// Check if this cluster exists
		if (heVList.find(std::make_pair(heIndex, v)) == heVList.end())
			continue;

		// Compute the distances
		heDistance = getHeDistance(heIndex);
		vDistance = getVDistance(v);

		// Add the concentration of each cluster in the group times its number of helium
		conc += getConcentration(heDistance, vDistance);
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
	setReactionConnectivity(heMomId);
	setDissociationConnectivity(heMomId);
	setReactionConnectivity(vMomId);
	setDissociationConnectivity(vMomId);

	// Visit all the reacting pairs
	std::for_each(effReactingList.begin(), effReactingList.end(),
			[this](ProductionPairMap::value_type const& currMapItem) {
				// The cluster is connecting to both clusters in the pair
				auto const& currPair = currMapItem.second;
				setReactionConnectivity(currPair.first.getId());
				setReactionConnectivity(currPair.first.getHeMomentumId());
				setReactionConnectivity(currPair.first.getVMomentumId());
				setReactionConnectivity(currPair.second.getId());
				setReactionConnectivity(currPair.second.getHeMomentumId());
				setReactionConnectivity(currPair.second.getVMomentumId());
			});

	// Visit all the combining pairs
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[this](CombiningClusterMap::value_type const& currMapItem) {
				// The cluster is connecting to the combining cluster
				auto const& currComb = currMapItem.second;
				setReactionConnectivity(currComb.first.getId());
				setReactionConnectivity(currComb.first.getHeMomentumId());
				setReactionConnectivity(currComb.first.getVMomentumId());
			});

	// Loop over all the dissociating pairs
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this](DissociationPairMap::value_type const& currMapItem) {
				// The cluster is connecting to the combining cluster
				auto const& currPair = currMapItem.second;
				setDissociationConnectivity(currPair.first.getId());
				setDissociationConnectivity(currPair.first.getHeMomentumId());
				setDissociationConnectivity(currPair.first.getVMomentumId());
			});

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	// Initialize the partial vector for the momentum
	int dof = network.getDOF();
	psiHeMomentPartials.resize(dof, 0.0);
	psiVMomentPartials.resize(dof, 0.0);

	return;
}

double PSISuperCluster::getDissociationFlux() {
	// Initial declarations
	double flux = 0.0;

	// Sum over all the dissociating pairs
	// TODO consider using std::accumulate.  May also want to change side
	// effect of updating member variables heMomentumFlux and
	// vMomentumFlux here.
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this,&flux](DissociationPairMap::value_type const& currMapItem) {
				auto const& currPair = currMapItem.second;

				// Get the dissociating clusters
				auto const& dissociatingCluster = currPair.first;
				double l0A = dissociatingCluster.getConcentration(0.0, 0.0);
				double lHeA = dissociatingCluster.getHeMomentum();
				double lVA = dissociatingCluster.getVMomentum();
				// Update the flux
				auto value = currPair.kConstant / (double) nTot;
				flux += value * (currPair.a00 * l0A + currPair.a10 * lHeA + currPair.a20 * lVA);
				// Compute the momentum fluxes
				heMomentumFlux += value
				* (currPair.a01 * l0A + currPair.a11 * lHeA + currPair.a21 * lVA);
				vMomentumFlux += value
				* (currPair.a02 * l0A + currPair.a12 * lHeA + currPair.a22 * lVA);
			});

	// Return the flux
	return flux;
}

double PSISuperCluster::getEmissionFlux() {
	// Initial declarations
	double flux = 0.0;

	// Loop over all the emission pairs
	// TODO consider using std::accumulate.  May also want to change side
	// effect of updating member variables heMomentumFlux and
	// vMomentumFlux here.
	std::for_each(effEmissionList.begin(), effEmissionList.end(),
			[this,&flux](DissociationPairMap::value_type const& currMapItem) {
				auto const& currPair = currMapItem.second;

				// Update the flux
				auto value = currPair.kConstant / (double) nTot;
				flux += value * (currPair.a00 * l0 + currPair.a10 * l1He + currPair.a20 * l1V);
				// Compute the momentum fluxes
				heMomentumFlux -= value
				* (currPair.a01 * l0 + currPair.a11 * l1He + currPair.a21 * l1V);
				vMomentumFlux -= value
				* (currPair.a02 * l0 + currPair.a12 * l1He + currPair.a22 * l1V);
			});

	return flux;
}

double PSISuperCluster::getProductionFlux() {
	// Local declarations
	double flux = 0.0;

	// Sum over all the reacting pairs
	// TODO consider using std::accumulate.  May also want to change side
	// effect of updating member variables heMomentumFlux and
	// vMomentumFlux here.
	std::for_each(effReactingList.begin(), effReactingList.end(),
			[this,&flux](ProductionPairMap::value_type const& currMapItem) {

				auto const& currPair = currMapItem.second;

				// Get the two reacting clusters
				auto const& firstReactant = currPair.first;
				auto const& secondReactant = currPair.second;
				double l0A = firstReactant.getConcentration(0.0, 0.0);
				double l0B = secondReactant.getConcentration(0.0, 0.0);
				double lHeA = firstReactant.getHeMomentum();
				double lHeB = secondReactant.getHeMomentum();
				double lVA = firstReactant.getVMomentum();
				double lVB = secondReactant.getVMomentum();
				// Update the flux
				auto value = currPair.kConstant / (double) nTot;
				flux += value
				* (currPair.a000 * l0A * l0B + currPair.a010 * l0A * lHeB
						+ currPair.a020 * l0A * lVB + currPair.a100 * lHeA * l0B
						+ currPair.a110 * lHeA * lHeB + currPair.a120 * lHeA * lVB
						+ currPair.a200 * lVA * l0B + currPair.a210 * lVA * lHeB
						+ currPair.a220 * lVA * lVB);
				// Compute the momentum fluxes
				heMomentumFlux += value
				* (currPair.a001 * l0A * l0B + currPair.a011 * l0A * lHeB
						+ currPair.a021 * l0A * lVB + currPair.a101 * lHeA * l0B
						+ currPair.a111 * lHeA * lHeB + currPair.a121 * lHeA * lVB
						+ currPair.a201 * lVA * l0B + currPair.a211 * lVA * lHeB
						+ currPair.a221 * lVA * lVB);
				vMomentumFlux += value
				* (currPair.a002 * l0A * l0B + currPair.a012 * l0A * lHeB
						+ currPair.a022 * l0A * lVB + currPair.a102 * lHeA * l0B
						+ currPair.a112 * lHeA * lHeB + currPair.a122 * lHeA * lVB
						+ currPair.a202 * lVA * l0B + currPair.a212 * lVA * lHeB
						+ currPair.a222 * lVA * lVB);
			});

	// Return the production flux
	return flux;
}

double PSISuperCluster::getCombinationFlux() {
	// Local declarations
	double flux = 0.0;

	// Sum over all the combining clusters
	// TODO consider using std::accumulate.  May also want to change side
	// effect of updating member variables heMomentumFlux and
	// vMomentumFlux here.
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[this,&flux](CombiningClusterMap::value_type const& currMapItem) {
				// Get the combining cluster
				auto const& currComb = currMapItem.second;
				auto const& combiningCluster = currComb.first;
				double l0B = combiningCluster.getConcentration(0.0, 0.0);
				double lHeB = combiningCluster.getHeMomentum();
				double lVB = combiningCluster.getVMomentum();
				// Update the flux
				auto value = currComb.kConstant / (double) nTot;
				flux += value
				* (currComb.a000 * l0B * l0 + currComb.a100 * l0B * l1He
						+ currComb.a200 * l0B * l1V + currComb.a010 * lHeB * l0
						+ currComb.a110 * lHeB * l1He + currComb.a210 * lHeB * l1V
						+ currComb.a020 * lVB * l0 + currComb.a120 * lVB * l1He
						+ currComb.a220 * lVB * l1V);
				// Compute the momentum fluxes
				heMomentumFlux -= value
				* (currComb.a001 * l0B * l0 + currComb.a101 * l0B * l1He
						+ currComb.a201 * l0B * l1V + currComb.a011 * lHeB * l0
						+ currComb.a111 * lHeB * l1He + currComb.a211 * lHeB * l1V
						+ currComb.a021 * lVB * l0 + currComb.a121 * lVB * l1He
						+ currComb.a221 * lVB * l1V);
				vMomentumFlux -= value
				* (currComb.a002 * l0B * l0 + currComb.a102 * l0B * l1He
						+ currComb.a202 * l0B * l1V + currComb.a012 * lHeB * l0
						+ currComb.a112 * lHeB * l1He + currComb.a212 * lHeB * l1V
						+ currComb.a022 * lVB * l0 + currComb.a122 * lVB * l1He
						+ currComb.a222 * lVB * l1V);
			});

	return flux;
}

void PSISuperCluster::getPartialDerivatives(
		std::vector<double> & partials) const {
	// Reinitialize the momentum partial derivatives vector
	std::fill(psiHeMomentPartials.begin(), psiHeMomentPartials.end(), 0.0);
	std::fill(psiVMomentPartials.begin(), psiVMomentPartials.end(), 0.0);

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
				double l0A = firstReactant.getConcentration(0.0, 0.0);
				double l0B = secondReactant.getConcentration(0.0, 0.0);
				double lHeA = firstReactant.getHeMomentum();
				double lHeB = secondReactant.getHeMomentum();
				double lVA = firstReactant.getVMomentum();
				double lVB = secondReactant.getVMomentum();

				// Compute the contribution from the first part of the reacting pair
				auto value = currPair.kConstant / (double) nTot;
				auto index = firstReactant.getId() - 1;
				partials[index] += value
				* (currPair.a000 * l0B + currPair.a010 * lHeB + currPair.a020 * lVB);
				psiHeMomentPartials[index] += value
				* (currPair.a001 * l0B + currPair.a011 * lHeB + currPair.a021 * lVB);
				psiVMomentPartials[index] += value
				* (currPair.a002 * l0B + currPair.a012 * lHeB + currPair.a022 * lVB);
				index = firstReactant.getHeMomentumId() - 1;
				partials[index] += value
				* (currPair.a100 * l0B + currPair.a110 * lHeB + currPair.a120 * lVB);
				psiHeMomentPartials[index] += value
				* (currPair.a101 * l0B + currPair.a111 * lHeB + currPair.a121 * lVB);
				psiVMomentPartials[index] += value
				* (currPair.a102 * l0B + currPair.a112 * lHeB + currPair.a122 * lVB);
				index = firstReactant.getVMomentumId() - 1;
				partials[index] += value
				* (currPair.a200 * l0B + currPair.a210 * lHeB + currPair.a220 * lVB);
				psiHeMomentPartials[index] += value
				* (currPair.a201 * l0B + currPair.a211 * lHeB + currPair.a221 * lVB);
				psiVMomentPartials[index] += value
				* (currPair.a202 * l0B + currPair.a212 * lHeB + currPair.a222 * lVB);
				// Compute the contribution from the second part of the reacting pair
				index = secondReactant.getId() - 1;
				partials[index] += value
				* (currPair.a000 * l0A + currPair.a100 * lHeA + currPair.a200 * lVA);
				psiHeMomentPartials[index] += value
				* (currPair.a001 * l0A + currPair.a101 * lHeA + currPair.a201 * lVA);
				psiVMomentPartials[index] += value
				* (currPair.a002 * l0A + currPair.a102 * lHeA + currPair.a202 * lVA);
				index = secondReactant.getHeMomentumId() - 1;
				partials[index] += value
				* (currPair.a010 * l0A + currPair.a110 * lHeA + currPair.a210 * lVA);
				psiHeMomentPartials[index] += value
				* (currPair.a011 * l0A + currPair.a111 * lHeA + currPair.a211 * lVA);
				psiVMomentPartials[index] += value
				* (currPair.a012 * l0A + currPair.a112 * lHeA + currPair.a212 * lVA);
				index = secondReactant.getVMomentumId() - 1;
				partials[index] += value
				* (currPair.a020 * l0A + currPair.a120 * lHeA + currPair.a220 * lVA);
				psiHeMomentPartials[index] += value
				* (currPair.a021 * l0A + currPair.a121 * lHeA + currPair.a221 * lVA);
				psiVMomentPartials[index] += value
				* (currPair.a022 * l0A + currPair.a122 * lHeA + currPair.a222 * lVA);
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
				double l0B = cluster.getConcentration(0.0, 0.0);
				double lHeB = cluster.getHeMomentum();
				double lVB = cluster.getVMomentum();

				// Compute the contribution from the combining cluster
				auto value = currComb.kConstant / (double) nTot;
				auto index = cluster.getId() - 1;
				partials[index] -= value
				* (currComb.a000 * l0 + currComb.a100 * l1He + currComb.a200 * l1V);
				psiHeMomentPartials[index] -= value
				* (currComb.a001 * l0 + currComb.a101 * l1He + currComb.a201 * l1V);
				psiVMomentPartials[index] -= value
				* (currComb.a002 * l0 + currComb.a102 * l1He + currComb.a202 * l1V);
				index = cluster.getHeMomentumId() - 1;
				partials[index] -= value
				* (currComb.a010 * l0 + currComb.a110 * l1He + currComb.a210 * l1V);
				psiHeMomentPartials[index] -= value
				* (currComb.a011 * l0 + currComb.a111 * l1He + currComb.a211 * l1V);
				psiVMomentPartials[index] -= value
				* (currComb.a012 * l0 + currComb.a112 * l1He + currComb.a212 * l1V);
				index = cluster.getVMomentumId() - 1;
				partials[index] -= value
				* (currComb.a020 * l0 + currComb.a120 * l1He + currComb.a220 * l1V);
				psiHeMomentPartials[index] -= value
				* (currComb.a021 * l0 + currComb.a121 * l1He + currComb.a221 * l1V);
				psiVMomentPartials[index] -= value
				* (currComb.a022 * l0 + currComb.a122 * l1He + currComb.a222 * l1V);
				// Compute the contribution from this cluster
				index = id - 1;
				partials[index] -= value
				* (currComb.a000 * l0B + currComb.a010 * lHeB + currComb.a020 * lVB);
				psiHeMomentPartials[index] -= value
				* (currComb.a001 * l0B + currComb.a011 * lHeB + currComb.a021 * lVB);
				psiVMomentPartials[index] -= value
				* (currComb.a002 * l0B + currComb.a012 * lHeB + currComb.a022 * lVB);
				index = heMomId - 1;
				partials[index] -= value
				* (currComb.a100 * l0B + currComb.a110 * lHeB + currComb.a120 * lVB);
				psiHeMomentPartials[index] -= value
				* (currComb.a101 * l0B + currComb.a111 * lHeB + currComb.a121 * lVB);
				psiVMomentPartials[index] -= value
				* (currComb.a102 * l0B + currComb.a112 * lHeB + currComb.a122 * lVB);
				index = vMomId - 1;
				partials[index] -= value
				* (currComb.a200 * l0B + currComb.a210 * lHeB + currComb.a220 * lVB);
				psiHeMomentPartials[index] -= value
				* (currComb.a201 * l0B + currComb.a211 * lHeB + currComb.a221 * lVB);
				psiVMomentPartials[index] -= value
				* (currComb.a202 * l0B + currComb.a212 * lHeB + currComb.a222 * lVB);
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
				auto index = cluster.getId() - 1;
				partials[index] += value * (currPair.a00);
				psiHeMomentPartials[index] += value * (currPair.a01);
				psiVMomentPartials[index] += value * (currPair.a02);
				index = cluster.getHeMomentumId() - 1;
				partials[index] += value * (currPair.a10);
				psiHeMomentPartials[index] += value * (currPair.a11);
				psiVMomentPartials[index] += value * (currPair.a12);
				index = cluster.getVMomentumId() - 1;
				partials[index] += value * (currPair.a20);
				psiHeMomentPartials[index] += value * (currPair.a21);
				psiVMomentPartials[index] += value * (currPair.a22);
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
				auto index = id - 1;
				partials[index] -= value * (currPair.a00);
				psiHeMomentPartials[index] -= value * (currPair.a01);
				psiVMomentPartials[index] -= value * (currPair.a02);
				index = heMomId - 1;
				partials[index] -= value * (currPair.a10);
				psiHeMomentPartials[index] -= value * (currPair.a11);
				psiVMomentPartials[index] -= value * (currPair.a12);
				index = vMomId - 1;
				partials[index] -= value * (currPair.a20);
				psiHeMomentPartials[index] -= value * (currPair.a21);
				psiVMomentPartials[index] -= value * (currPair.a22);
			});

	return;
}

void PSISuperCluster::getHeMomentPartialDerivatives(
		std::vector<double> & partials) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = psiHeMomentPartials[i];
	}

	return;
}

void PSISuperCluster::getVMomentPartialDerivatives(
		std::vector<double> & partials) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = psiVMomentPartials[i];
	}

	return;
}

void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::ProductionCoefficientBase const& curr) const {

	os << "a[0-2][0-2][0-2]:" << ' ' << curr.a000 << ' ' << curr.a001 << ' '
			<< curr.a002 << ' ' << curr.a100 << ' ' << curr.a101 << ' '
			<< curr.a102 << ' ' << curr.a200 << ' ' << curr.a201 << ' '
			<< curr.a202 << ' ' << curr.a010 << ' ' << curr.a011 << ' '
			<< curr.a012 << ' ' << curr.a020 << ' ' << curr.a021 << ' '
			<< curr.a022 << ' ' << curr.a110 << ' ' << curr.a111 << ' '
			<< curr.a112 << ' ' << curr.a120 << ' ' << curr.a121 << ' '
			<< curr.a122 << ' ' << curr.a210 << ' ' << curr.a211 << ' '
			<< curr.a212 << ' ' << curr.a220 << ' ' << curr.a221 << ' '
			<< curr.a222;
}

void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::SuperClusterDissociationPair const& currPair) const {

	os << "a[0-2][0-2]:" << ' ' << currPair.a00 << ' ' << currPair.a01 << ' '
			<< currPair.a02 << ' ' << currPair.a10 << ' ' << currPair.a11 << ' '
			<< currPair.a12 << ' ' << currPair.a20 << ' ' << currPair.a21 << ' '
			<< currPair.a22;
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

