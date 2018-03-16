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

PSISuperCluster::PSISuperCluster(double numHe, double numV, int _nTot,
		int heWidth, int vWidth, IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(_network, registry, buildName(numHe, numV)), nTot(_nTot), heBounds(
				0, 0), vBounds(0, 0), l0(0.0), l1He(0.0), l1V(0.0), heMomentumFlux(
				0.0), vMomentumFlux(0.0) {
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	numAtom[0] = numHe;
	numAtom[1] = numV;
	size = (int) (numHe + numV);

	// Update the composition map
	composition[toCompIdx(Species::He)] = (int) numHe;
	composition[toCompIdx(Species::V)] = (int) numV;

	// Set the width
	sectionWidth[0] = heWidth;
	sectionWidth[1] = vWidth;

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
	double firstDistance[3] = { 1.0, 0.0, 0.0 }, secondDistance[3] = { 1.0, 0.0,
			0.0 };
	if (prodPair.first.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(prodPair.first);
		firstDistance[1] = super.getDistance(c, 0);
		firstDistance[2] = super.getDistance(d, 1);
	}
	if (prodPair.second.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(prodPair.second);
		secondDistance[1] = super.getDistance(c, 0);
		secondDistance[2] = super.getDistance(d, 1);
	}
	double factor[3] = { 1.0, 0.0, 0.0 };
	factor[1] = getFactor(a, 0);
	factor[2] = getFactor(b, 1);
	// First is A, second is B, in A + B -> this
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < 3; j++) {
			for (int i = 0; i < 3; i++) {
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

				// Use names corresponding to those in single version.
				int a = currPRI.numHe;
				int b = currPRI.numV;
				int c = currPRI.i;
				int d = currPRI.j;

				// Update the coefficients
				double firstDistance[3] = {1.0, 0.0, 0.0}, secondDistance[3] = {1.0, 0.0, 0.0};
				if (prodPair.first.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster const&>(prodPair.first);
					firstDistance[1] = super.getDistance(c, 0);
					firstDistance[2] = super.getDistance(d, 1);
				}
				if (prodPair.second.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster const&>(prodPair.second);
					secondDistance[1] = super.getDistance(c, 0);
					secondDistance[2] = super.getDistance(d, 1);
				}
				double factor[3] = {1.0, 0.0, 0.0};
				factor[1] = getFactor(a, 0);
				factor[2] = getFactor(b, 1);
				// First is A, second is B, in A + B -> this
				for (int k = 0; k < 3; k++) {
					for (int j = 0; j < 3; j++) {
						for (int i = 0; i < 3; i++) {
							prodPair.coefs[i][j][k] += firstDistance[i] * secondDistance[j] * factor[k];
						}
					}
				}
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
	double distance[3] = { 1.0, 0.0, 0.0 };
	double factor[3] = { 1.0, 0.0, 0.0 };
	distance[1] = getDistance(a, 0);
	distance[2] = getDistance(b, 1);
	factor[1] = getFactor(a, 0);
	factor[2] = getFactor(b, 1);
	// This is A, itBis is B, in A + B -> C
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < 3; j++) {
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
			[this,&combCluster](const PendingProductionReactionInfo& currPRInfo) {

				// Use names corresponding to the single-item version.
				int a = currPRInfo.i;
				int b = currPRInfo.j;

				// Update the coefficients
				double distance[3] = {1.0, 0.0, 0.0};
				double factor[3] = {1.0, 0.0, 0.0};
				distance[1] = getDistance(a, 0);
				distance[2] = getDistance(b, 1);
				factor[1] = getFactor(a, 0);
				factor[2] = getFactor(b, 1);
				// This is A, itBis is B, in A + B -> C
				for (int k = 0; k < 3; k++) {
					for (int j = 0; j < 3; j++) {
						combCluster.coefs[j][0][k] += distance[j] * factor[k];
					}
				}
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
	double distance[3] = { 1.0, 0.0, 0.0 };
	if (reaction.dissociating.getType() == ReactantType::PSISuper) {
		auto const& super =
				static_cast<PSICluster const&>(reaction.dissociating);
		distance[1] = super.getDistance(a, 0);
		distance[2] = super.getDistance(b, 1);
	}
	double factor[3] = { 1.0, 0.0, 0.0 };
	factor[1] = getFactor(c, 0);
	factor[2] = getFactor(d, 1);

	// A is the dissociating cluster
	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < 3; i++) {
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

				// Use names corresponding to the single-item version.
				int a = currPRI.numHe;
				int b = currPRI.numV;
				int c = currPRI.i;
				int d = currPRI.j;

				// Update the coefficients
				double distance[3] = {1.0, 0.0, 0.0};
				if (reaction.dissociating.getType() == ReactantType::PSISuper) {
					auto const& super =
					static_cast<PSICluster const&>(reaction.dissociating);
					distance[1] = super.getDistance(a, 0);
					distance[2] = super.getDistance(b, 1);
				}
				double factor[3] = {1.0, 0.0, 0.0};
				factor[1] = getFactor(c, 0);
				factor[2] = getFactor(d, 1);

				// A is the dissociating cluster
				for (int j = 0; j < 3; j++) {
					for (int i = 0; i < 3; i++) {
						dissPair.coefs[i][j] += distance[i] * factor[j];
					}
				}
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

	// Update the coefficients
	double distance[3] = { 1.0, 0.0, 0.0 };
	double factor[3] = { 1.0, 0.0, 0.0 };
	distance[1] = getDistance(a, 0);
	distance[2] = getDistance(b, 1);
	factor[1] = getFactor(a, 0);
	factor[2] = getFactor(b, 1);
	// A is the dissociating cluster
	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < 3; i++) {
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

				// Use same names as used in single version.
				int a = currPRI.numHe;
				int b = currPRI.numV;

				// Update the coefficients
				double distance[3] = {1.0, 0.0, 0.0};
				double factor[3] = {1.0, 0.0, 0.0};
				distance[1] = getDistance(a, 0);
				distance[2] = getDistance(b, 1);
				factor[1] = getFactor(a, 0);
				factor[2] = getFactor(b, 1);
				// A is the dissociating cluster
				for (int j = 0; j < 3; j++) {
					for (int i = 0; i < 3; i++) {
						dissPair.coefs[i][j] += distance[i] * factor[j];
					}
				}
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
	if (sectionWidth[0] == 1)
		dispersion[0] = 1.0;
	else
		dispersion[0] = 2.0
				* (nHeSquare - (numAtom[0] * (double) nTot * numAtom[0]))
				/ ((double) (nTot * (sectionWidth[0] - 1)));

	if (sectionWidth[1] == 1)
		dispersion[1] = 1.0;
	else
		dispersion[1] = 2.0
				* (nVSquare - (numAtom[1] * (double) nTot * numAtom[1]))
				/ ((double) (nTot * (sectionWidth[1] - 1)));

	// Set the boundaries
	heBounds = IntegerRange<IReactant::SizeType>(
			static_cast<IReactant::SizeType>((numAtom[0]
					- (double) sectionWidth[0] / 2.0) + 1),
			static_cast<IReactant::SizeType>((numAtom[0]
					- (double) sectionWidth[0] / 2.0) + sectionWidth[0]) + 1);
	vBounds = IntegerRange<IReactant::SizeType>(
			static_cast<IReactant::SizeType>((numAtom[1]
					- (double) sectionWidth[1] / 2.0) + 1),
			static_cast<IReactant::SizeType>((numAtom[1]
					- (double) sectionWidth[1] / 2.0) + sectionWidth[1]) + 1);

	return;
}

double PSISuperCluster::getTotalConcentration() const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
	for (auto const& pair : heVList) {
		// Compute the distances
		heDistance = getDistance(pair.first, 0);
		vDistance = getDistance(pair.second, 1);

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
		heDistance = getDistance(pair.first, 0);
		vDistance = getDistance(pair.second, 1);

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
		heDistance = getDistance(pair.first, 0);
		vDistance = getDistance(pair.second, 1);

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
	for (int j = 0; j < sectionWidth[0]; j++) {
		// Compute the helium index
		heIndex = (int) (numAtom[0] - (double) sectionWidth[0] / 2.0) + j + 1;

		// Check if this cluster exists
		if (heVList.find(std::make_pair(heIndex, v)) == heVList.end())
			continue;

		// Compute the distances
		heDistance = getDistance(heIndex, 0);
		vDistance = getDistance(v, 1);

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
				double lA[3] = {};
				lA[0] = dissociatingCluster.getConcentration(0.0, 0.0);
				lA[1] = dissociatingCluster.getHeMomentum();
				lA[2] = dissociatingCluster.getVMomentum();

				double sum[3] = {};
				for (int j = 0; j < 3; j++) {
					for (int i = 0; i < 3; i++) {
						sum[j] += currPair.coefs[i][j] * lA[i];
					}
				}
				// Update the flux
				auto value = currPair.kConstant / (double) nTot;
				flux += value * sum[0];
				// Compute the momentum fluxes
				heMomentumFlux += value * sum[1];
				vMomentumFlux += value * sum[2];
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
				double lA[3] = {};
				lA[0] = l0;
				lA[1] = l1He;
				lA[2] = l1V;

				double sum[3] = {};
				for (int j = 0; j < 3; j++) {
					for (int i = 0; i < 3; i++) {
						sum[j] += currPair.coefs[i][j] * lA[i];
					}
				}
				// Update the flux
				auto value = currPair.kConstant / (double) nTot;
				flux += value * sum[0];
				// Compute the momentum fluxes
				heMomentumFlux -= value * sum[1];
				vMomentumFlux -= value * sum[2];
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
				double lA[3] = {}, lB[3] = {};
				lA[0] = firstReactant.getConcentration(0.0, 0.0);
				lB[0] = secondReactant.getConcentration(0.0, 0.0);
				lA[1] = firstReactant.getHeMomentum();
				lB[1] = secondReactant.getHeMomentum();
				lA[2] = firstReactant.getVMomentum();
				lB[2] = secondReactant.getVMomentum();

				double sum[3] = {};
				for (int k = 0; k < 3; k++) {
					for (int j = 0; j < 3; j++) {
						for (int i = 0; i < 3; i++) {
							sum[k] += currPair.coefs[j][i][k] * lA[j] * lB[i];
						}
					}
				}

				// Update the flux
				auto value = currPair.kConstant / (double) nTot;
				flux += value * sum[0];
				// Compute the momentum fluxes
				heMomentumFlux += value * sum[1];
				vMomentumFlux += value * sum[2];
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
				double lA[3] = {}, lB[3] = {};
				lA[0] = l0;
				lB[0] = combiningCluster.getConcentration(0.0, 0.0);
				lA[1] = l1He;
				lB[1] = combiningCluster.getHeMomentum();
				lA[2] = l1V;
				lB[2] = combiningCluster.getVMomentum();

				double sum[3] = {};
				for (int k = 0; k < 3; k++) {
					for (int j = 0; j < 3; j++) {
						for (int i = 0; i < 3; i++) {
							sum[k] += currComb.coefs[i][j][k] * lA[i] * lB[j];
						}
					}
				}
				double l0B = combiningCluster.getConcentration(0.0, 0.0);
				double lHeB = combiningCluster.getHeMomentum();
				double lVB = combiningCluster.getVMomentum();
				// Update the flux
				auto value = currComb.kConstant / (double) nTot;
				flux += value * sum[0];
				// Compute the momentum fluxes
				heMomentumFlux -= value * sum[1];
				vMomentumFlux -= value * sum[2];
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
				double lA[3] = {}, lB[3] = {};
				lA[0] = firstReactant.getConcentration(0.0, 0.0);
				lB[0] = secondReactant.getConcentration(0.0, 0.0);
				lA[1] = firstReactant.getHeMomentum();
				lB[1] = secondReactant.getHeMomentum();
				lA[2] = firstReactant.getVMomentum();
				lB[2] = secondReactant.getVMomentum();

				double sum[3][3][2] = {};
				for (int k = 0; k < 3; k++) {
					for (int j = 0; j < 3; j++) {
						for (int i = 0; i < 3; i++) {
							sum[k][j][0] += currPair.coefs[j][i][k] * lB[i];
							sum[k][j][1] += currPair.coefs[i][j][k] * lA[i];
						}
					}
				}

				// Compute the contribution from the first part of the reacting pair
				auto value = currPair.kConstant / (double) nTot;
				auto index = firstReactant.getId() - 1;
				partials[index] += value * sum[0][0][0];
				psiHeMomentPartials[index] += value * sum[1][0][0];
				psiVMomentPartials[index] += value * sum[2][0][0];
				index = firstReactant.getHeMomentumId() - 1;
				partials[index] += value * sum[0][1][0];
				psiHeMomentPartials[index] += value * sum[1][1][0];
				psiVMomentPartials[index] += value * sum[2][1][0];
				index = firstReactant.getVMomentumId() - 1;
				partials[index] += value * sum[0][2][0];
				psiHeMomentPartials[index] += value * sum[1][2][0];
				psiVMomentPartials[index] += value * sum[2][2][0];
				// Compute the contribution from the second part of the reacting pair
				index = secondReactant.getId() - 1;
				partials[index] += value * sum[0][0][1];
				psiHeMomentPartials[index] += value * sum[1][0][1];
				psiVMomentPartials[index] += value * sum[2][0][1];
				index = secondReactant.getHeMomentumId() - 1;
				partials[index] += value * sum[0][1][1];
				psiHeMomentPartials[index] += value * sum[1][1][1];
				psiVMomentPartials[index] += value * sum[2][1][1];
				index = secondReactant.getVMomentumId() - 1;
				partials[index] += value * sum[0][2][1];
				psiHeMomentPartials[index] += value * sum[1][2][1];
				psiVMomentPartials[index] += value * sum[2][2][1];
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
				double lA[3] = {}, lB[3] = {};
				lA[0] = l0;
				lB[0] = cluster.getConcentration(0.0, 0.0);
				lA[1] = l1He;
				lB[1] = cluster.getHeMomentum();
				lA[2] = l1V;
				lB[2] = cluster.getVMomentum();

				double sum[3][3][2] = {};
				for (int k = 0; k < 3; k++) {
					for (int j = 0; j < 3; j++) {
						for (int i = 0; i < 3; i++) {
							sum[k][j][0] += currComb.coefs[i][j][k] * lA[i];
							sum[k][j][1] += currComb.coefs[j][i][k] * lB[i];
						}
					}
				}

				// Compute the contribution from the combining cluster
				auto value = currComb.kConstant / (double) nTot;
				auto index = cluster.getId() - 1;
				partials[index] -= value * sum[0][0][0];
				psiHeMomentPartials[index] -= value * sum[1][0][0];
				psiVMomentPartials[index] -= value * sum[2][0][0];
				index = cluster.getHeMomentumId() - 1;
				partials[index] -= value * sum[0][1][0];
				psiHeMomentPartials[index] -= value * sum[1][1][0];
				psiVMomentPartials[index] -= value * sum[2][1][0];
				index = cluster.getVMomentumId() - 1;
				partials[index] -= value * sum[0][2][0];
				psiHeMomentPartials[index] -= value * sum[1][2][0];
				psiVMomentPartials[index] -= value * sum[2][2][0];
				// Compute the contribution from this cluster
				index = id - 1;
				partials[index] -= value * sum[0][0][1];
				psiHeMomentPartials[index] -= value * sum[1][0][1];
				psiVMomentPartials[index] -= value * sum[2][0][1];
				index = heMomId - 1;
				partials[index] -= value * sum[0][1][1];
				psiHeMomentPartials[index] -= value * sum[1][1][1];
				psiVMomentPartials[index] -= value * sum[2][1][1];
				index = vMomId - 1;
				partials[index] -= value * sum[0][2][1];
				psiHeMomentPartials[index] -= value * sum[1][2][1];
				psiVMomentPartials[index] -= value * sum[2][2][1];
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
				partials[index] += value * (currPair.coefs[0][0]);
				psiHeMomentPartials[index] += value * (currPair.coefs[0][1]);
				psiVMomentPartials[index] += value * (currPair.coefs[0][2]);
				index = cluster.getHeMomentumId() - 1;
				partials[index] += value * (currPair.coefs[1][0]);
				psiHeMomentPartials[index] += value * (currPair.coefs[1][1]);
				psiVMomentPartials[index] += value * (currPair.coefs[1][2]);
				index = cluster.getVMomentumId() - 1;
				partials[index] += value * (currPair.coefs[2][0]);
				psiHeMomentPartials[index] += value * (currPair.coefs[2][1]);
				psiVMomentPartials[index] += value * (currPair.coefs[2][2]);
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
				partials[index] -= value * (currPair.coefs[0][0]);
				psiHeMomentPartials[index] -= value * (currPair.coefs[0][1]);
				psiVMomentPartials[index] -= value * (currPair.coefs[0][2]);
				index = heMomId - 1;
				partials[index] -= value * (currPair.coefs[1][0]);
				psiHeMomentPartials[index] -= value * (currPair.coefs[1][1]);
				psiVMomentPartials[index] -= value * (currPair.coefs[1][2]);
				index = vMomId - 1;
				partials[index] -= value * (currPair.coefs[2][0]);
				psiHeMomentPartials[index] -= value * (currPair.coefs[2][1]);
				psiVMomentPartials[index] -= value * (currPair.coefs[2][2]);
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

	os << "a[0-2][0-2][0-2]:" << ' ' << curr.coefs[0][0][0] << ' '
			<< curr.coefs[0][0][1] << ' ' << curr.coefs[0][0][2] << ' '
			<< curr.coefs[1][0][0] << ' ' << curr.coefs[1][0][1] << ' '
			<< curr.coefs[1][0][2] << ' ' << curr.coefs[2][0][0] << ' '
			<< curr.coefs[2][0][1] << ' ' << curr.coefs[2][0][2] << ' '
			<< curr.coefs[0][1][0] << ' ' << curr.coefs[0][1][1] << ' '
			<< curr.coefs[0][1][2] << ' ' << curr.coefs[0][2][0] << ' '
			<< curr.coefs[0][2][1] << ' ' << curr.coefs[0][2][2] << ' '
			<< curr.coefs[1][1][0] << ' ' << curr.coefs[1][1][1] << ' '
			<< curr.coefs[1][1][2] << ' ' << curr.coefs[1][2][0] << ' '
			<< curr.coefs[1][2][1] << ' ' << curr.coefs[1][2][2] << ' '
			<< curr.coefs[2][1][0] << ' ' << curr.coefs[2][1][1] << ' '
			<< curr.coefs[2][1][2] << ' ' << curr.coefs[2][2][0] << ' '
			<< curr.coefs[2][2][1] << ' ' << curr.coefs[2][2][2];
}

void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::SuperClusterDissociationPair const& currPair) const {

	os << "a[0-2][0-2]:" << ' ' << currPair.coefs[0][0] << ' '
			<< currPair.coefs[0][1] << ' ' << currPair.coefs[0][2] << ' '
			<< currPair.coefs[1][0] << ' ' << currPair.coefs[1][1] << ' '
			<< currPair.coefs[1][2] << ' ' << currPair.coefs[2][0] << ' '
			<< currPair.coefs[2][1] << ' ' << currPair.coefs[2][2];
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

