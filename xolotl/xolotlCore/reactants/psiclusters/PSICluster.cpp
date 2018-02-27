#include <algorithm>
#include <cassert>
#include "PSICluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

void PSICluster::resultFrom(ProductionReaction& reaction, int a, int b, int c,
		int d) {

	// Add a cluster pair for the given reaction.
	reactingPairs.emplace_back(reaction,
			static_cast<PSICluster&>(reaction.first),
			static_cast<PSICluster&>(reaction.second));
	auto& newPair = reactingPairs.back();

	// NB: newPair's reactants are same as reaction's.
	// So use newPair only from here on.
	// TODO Any way to enforce this beyond splitting it into two functions?

	// Update the coefficients
	double firstHeDistance = 0.0, firstVDistance = 0.0, secondHeDistance = 0.0,
			secondVDistance = 0.0;
	if (newPair.first.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(newPair.first);
		firstHeDistance = super.getHeDistance(c);
		firstVDistance = super.getVDistance(d);
	}
	if (newPair.second.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(newPair.second);
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

	return;
}

void PSICluster::resultFrom(ProductionReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& prInfos) {

	// Add a cluster pair for the given reaction.
	reactingPairs.emplace_back(reaction,
			static_cast<PSICluster&>(reaction.first),
			static_cast<PSICluster&>(reaction.second));
	auto& newPair = reactingPairs.back();

	// NB: newPair's reactants are same as reaction's.
	// So use newPair only from here on.
	// TODO Any way to enforce this beyond splitting it into two functions?

	// Update the coefficients
	std::for_each(prInfos.begin(), prInfos.end(),
			[&newPair](const PendingProductionReactionInfo& currPRI) {

				// Use names that correspond to single version.
				int a = currPRI.numHe;
				int b = currPRI.numV;
				int c = currPRI.i;
				int d = currPRI.j;

				double firstHeDistance = 0.0, firstVDistance = 0.0, secondHeDistance = 0.0,
				secondVDistance = 0.0;
				if (newPair.first.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster const&>(newPair.first);
					firstHeDistance = super.getHeDistance(c);
					firstVDistance = super.getVDistance(d);
				}
				if (newPair.second.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster const&>(newPair.second);
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

void PSICluster::participateIn(ProductionReaction& reaction, int a, int b) {
	// Look for the other cluster
	auto& otherCluster = static_cast<PSICluster&>(
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
	if (otherCluster.getType() == ReactantType::PSISuper) {
		heDistance = otherCluster.getHeDistance(a);
		vDistance = otherCluster.getVDistance(b);
	}
	(*it).a0 += 1.0;
	(*it).a1 += heDistance;
	(*it).a2 += vDistance;

	return;
}

void PSICluster::participateIn(ProductionReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& prInfos) {
	// Look for the other cluster
	auto& otherCluster = static_cast<PSICluster&>(
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
			[&otherCluster,&combCluster](const PendingProductionReactionInfo& currPRInfo) {

				// Use names corresponding to the single-item version.
				int a = currPRInfo.i;
				int b = currPRInfo.j;

				double heDistance = 0.0, vDistance = 0.0;
				if (otherCluster.getType() == ReactantType::PSISuper) {
					heDistance = otherCluster.getHeDistance(a);
					vDistance = otherCluster.getVDistance(b);
				}
				combCluster.a0 += 1.0;
				combCluster.a1 += heDistance;
				combCluster.a2 += vDistance;
			});

	return;
}

void PSICluster::participateIn(DissociationReaction& reaction, int a, int b,
		int c, int d) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<PSICluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if the reaction was already added
	auto it =
			std::find_if(dissociatingPairs.rbegin(), dissociatingPairs.rend(),
					[&reaction,&emittedCluster](const ClusterPair& currPair) {
						return (&(reaction.dissociating) == &static_cast<PSICluster&>(currPair.first)) and
						(&emittedCluster == &static_cast<PSICluster&>(currPair.second));

					});
	if (it == dissociatingPairs.rend()) {

		// We did not already know about it.
		// Add the pair of them where it is important that the
		// dissociating cluster is the first one
		dissociatingPairs.emplace_back(reaction,
				static_cast<PSICluster&>(reaction.dissociating),
				static_cast<PSICluster&>(emittedCluster));
		it = dissociatingPairs.rbegin();
	}

	// Update the coefficients
	double firstHeDistance = 0.0, firstVDistance = 0.0;
	if (reaction.dissociating.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster&>(reaction.dissociating);
		firstHeDistance = super.getHeDistance(a);
		firstVDistance = super.getVDistance(b);
	}
	(*it).a00 += 1.0;
	(*it).a10 += firstHeDistance;
	(*it).a20 += firstVDistance;

	return;
}

void PSICluster::participateIn(DissociationReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& prInfos) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<PSICluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if the reaction was already added
	auto it =
			std::find_if(dissociatingPairs.rbegin(), dissociatingPairs.rend(),
					[&reaction,&emittedCluster](const ClusterPair& currPair) {
						return (&(reaction.dissociating) == &static_cast<PSICluster&>(currPair.first)) and
						(&emittedCluster == &static_cast<PSICluster&>(currPair.second));

					});
	if (it == dissociatingPairs.rend()) {

		// We did not already know about it.
		// Add the pair of them where it is important that the
		// dissociating cluster is the first one
		dissociatingPairs.emplace_back(reaction,
				static_cast<PSICluster&>(reaction.dissociating),
				static_cast<PSICluster&>(emittedCluster));
		it = dissociatingPairs.rbegin();
	}
	assert(it != dissociatingPairs.rend());
	auto& currPair = *it;

	// Update the coefficients
	std::for_each(prInfos.begin(), prInfos.end(),
			[&currPair,&reaction](const PendingProductionReactionInfo& currPRI) {

				// Use names corresponding to the single-item version.
				int a = currPRI.numHe;
				int b = currPRI.numV;

				double firstHeDistance = 0.0, firstVDistance = 0.0;
				if (reaction.dissociating.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster&>(reaction.dissociating);
					firstHeDistance = super.getHeDistance(a);
					firstVDistance = super.getVDistance(b);
				}
				currPair.a00 += 1.0;
				currPair.a10 += firstHeDistance;
				currPair.a20 += firstVDistance;
			});

	return;
}

void PSICluster::emitFrom(DissociationReaction& reaction, int a, int b, int c,
		int d) {

	// Note that we emit from the reaction's reactants according to
	// the given reaction.
	// TODO do we need to check to see whether we already know about
	// this reaction?
	emissionPairs.emplace_back(reaction,
			static_cast<PSICluster&>(reaction.first),
			static_cast<PSICluster&>(reaction.second));

	return;
}

void PSICluster::emitFrom(DissociationReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& prInfos) {

	// Note that we emit from the reaction's reactants according to
	// the given reaction.
	// TODO do we need to check to see whether we already know about
	// this reaction?
	emissionPairs.emplace_back(reaction,
			static_cast<PSICluster&>(reaction.first),
			static_cast<PSICluster&>(reaction.second));

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

std::vector<int> PSICluster::getReactionConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(reactionConnectivitySet, network.getDOF());
}

std::vector<int> PSICluster::getDissociationConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(dissociationConnectivitySet,
			network.getDOF());
}

void PSICluster::resetConnectivities() {
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
	setReactionConnectivity(heMomId);
	setDissociationConnectivity(heMomId);
	setReactionConnectivity(vMomId);
	setDissociationConnectivity(vMomId);

	// Loop on the effective reacting pairs
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[this](const ClusterPair& currPair) {
				// The cluster is connecting to both clusters in the pair
				setReactionConnectivity(currPair.first.id);
				setReactionConnectivity(currPair.second.id);
				setReactionConnectivity(currPair.first.heMomId);
				setReactionConnectivity(currPair.second.heMomId);
				setReactionConnectivity(currPair.first.vMomId);
				setReactionConnectivity(currPair.second.vMomId);
			});

	// Loop on the effective combining reactants
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[this](const CombiningCluster& cc) {
				// The cluster is connecting to the combining cluster
				setReactionConnectivity(cc.combining.id);
				setReactionConnectivity(cc.combining.heMomId);
				setReactionConnectivity(cc.combining.vMomId);
			});

	// Loop on the effective dissociating pairs
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[this](const ClusterPair& currPair) {
				// The cluster is connecting to the dissociating cluster which
				// is the first one by definition
				setDissociationConnectivity(currPair.first.id);
				setDissociationConnectivity(currPair.first.heMomId);
				setDissociationConnectivity(currPair.first.vMomId);
			});

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	return;
}

void PSICluster::updateFromNetwork() {

	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	return;
}

double PSICluster::getDissociationFlux() const {

	// Sum dissociation flux over all our dissociating clusters.
	double flux = std::accumulate(dissociatingPairs.begin(),
			dissociatingPairs.end(), 0.0,
			[](double running, const ClusterPair& currPair) {
				auto const& dissCluster = currPair.first;
				double l0A = dissCluster.getConcentration(0.0, 0.0);
				double lHeA = dissCluster.getHeMomentum();
				double lVA = dissCluster.getVMomentum();

				// Calculate the Dissociation flux
				return running +
				(currPair.reaction.kConstant *
						(currPair.a00 * l0A +
								currPair.a10 * lHeA +
								currPair.a20 * lVA));
			});

	// Return the flux
	return flux;
}

double PSICluster::getEmissionFlux() const {

	// Sum rate constants from all emission pair reactions.
	double flux = std::accumulate(emissionPairs.begin(), emissionPairs.end(),
			0.0, [](double running, const ClusterPair& currPair) {
				return running + currPair.reaction.kConstant;
			});

	return flux * concentration;
}

double PSICluster::getProductionFlux() const {

	// Sum production flux over all reacting pairs.
	double flux = std::accumulate(reactingPairs.begin(), reactingPairs.end(),
			0.0, [](double running, const ClusterPair& currPair) {

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
			return running + currPair.reaction.kConstant *
			(currPair.a00 * l0A * l0B + currPair.a01 * l0A * lHeB +
					currPair.a02 * l0A * lVB + currPair.a10 * lHeA * l0B +
					currPair.a11 * lHeA * lHeB + currPair.a12 * lHeA * lVB +
					currPair.a20 * lVA * l0B + currPair.a21 * lVA * lHeB +
					currPair.a22 * lVA * lVB);
		});

	// Return the production flux
	return flux;
}

double PSICluster::getCombinationFlux() const {

	// Sum combination flux over all clusters that combine with us.
	double flux = std::accumulate(combiningReactants.begin(),
			combiningReactants.end(), 0.0,
			[](double running, const CombiningCluster& cc) {

				// Get the cluster that combines with this one
				auto const& combiningCluster = cc.combining;
				double l0B = combiningCluster.getConcentration(0.0, 0.0);
				double lHeB = combiningCluster.getHeMomentum();
				double lVB = combiningCluster.getVMomentum();
				// Calculate the combination flux
				return running + (cc.reaction.kConstant *
						(cc.a0 * l0B + cc.a1 * lHeB + cc.a2 * lVB));

			});

	return flux * concentration;
}

std::vector<double> PSICluster::getPartialDerivatives() const {
	// Local Declarations
	std::vector<double> partials(network.getDOF(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return partials;
}

void PSICluster::getPartialDerivatives(std::vector<double> & partials) const {
	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return;
}

void PSICluster::getProductionPartialDerivatives(
		std::vector<double> & partials) const {

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[&partials](const ClusterPair& currPair) {
				// Get the two reacting clusters
				auto const& firstReactant = currPair.first;
				auto const& secondReactant = currPair.second;
				double l0A = firstReactant.getConcentration(0.0, 0.0);
				double l0B = secondReactant.getConcentration(0.0, 0.0);
				double lHeA = firstReactant.getHeMomentum();
				double lHeB = secondReactant.getHeMomentum();
				double lVA = firstReactant.getVMomentum();
				double lVB = secondReactant.getVMomentum();

				// Compute contribution from the first part of the reacting pair
				double value = currPair.reaction.kConstant;

				partials[firstReactant.id - 1] += value *
				(currPair.a00 * l0B + currPair.a01 * lHeB + currPair.a02 * lVB);

				partials[firstReactant.heMomId - 1] += value *
				(currPair.a10 * l0B + currPair.a11 * lHeB + currPair.a12 * lVB);

				partials[firstReactant.vMomId - 1] += value *
				(currPair.a20 * l0B + currPair.a21 * lHeB + currPair.a22 * lVB);

				// Compute contribution from the second part of the reacting pair
				partials[secondReactant.id - 1] += value *
				(currPair.a00 * l0A + currPair.a10 * lHeA + currPair.a20 * lVA);

				partials[secondReactant.heMomId - 1] += value *
				(currPair.a01 * l0A + currPair.a11 * lHeA + currPair.a21 * lVA);

				partials[secondReactant.vMomId - 1] += value *
				(currPair.a02 * l0A + currPair.a12 * lHeA + currPair.a22 * lVA);
			});

	return;
}

void PSICluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[this,&partials](const CombiningCluster& cc) {
				auto const& cluster = cc.combining;
				double l0B = cluster.getConcentration(0.0, 0.0);
				double lHeB = cluster.getHeMomentum();
				double lVB = cluster.getVMomentum();

				// Remember that the flux due to combinations is OUTGOING (-=)!
				// Compute the contribution from this cluster
				partials[id - 1] -= cc.reaction.kConstant
				* (cc.a0 * l0B + cc.a1 * lHeB + cc.a2 * lVB);
				// Compute the contribution from the combining cluster
				double value = cc.reaction.kConstant * concentration;
				partials[cluster.id - 1] -= value * cc.a0;
				partials[cluster.heMomId - 1] -= value * cc.a1;
				partials[cluster.vMomId - 1] -= value * cc.a2;
			});

	return;
}

void PSICluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[&partials](const ClusterPair& currPair) {
				// Get the dissociating cluster
				auto const& cluster = currPair.first;
				double value = currPair.reaction.kConstant;
				partials[cluster.id - 1] += value * currPair.a00;
				partials[cluster.heMomId - 1] += value * currPair.a10;
				partials[cluster.vMomId - 1] += value * currPair.a20;
			});

	return;
}

void PSICluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)
	double outgoingFlux = std::accumulate(emissionPairs.begin(),
			emissionPairs.end(), 0.0,
			[](double running, const ClusterPair& currPair) {
				return running + currPair.reaction.kConstant;
			});
	partials[id - 1] -= outgoingFlux;

	return;
}

void PSICluster::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

void PSICluster::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

double PSICluster::getLeftSideRate() const {

	// Sum rate constant-concentration product over combining reactants.
	double combiningRateTotal =
			std::accumulate(combiningReactants.begin(),
					combiningReactants.end(), 0.0,
					[](double running, const CombiningCluster& cc) {
						return running +
						(cc.reaction.kConstant * cc.combining.concentration * cc.a0);
					});

	// Sum rate constants over all emission pair reactions.
	double emissionRateTotal = std::accumulate(emissionPairs.begin(),
			emissionPairs.end(), 0.0,
			[](double running, const ClusterPair& currPair) {
				return running + (currPair.reaction.kConstant);
			});

	return combiningRateTotal + emissionRateTotal;
}

std::vector<int> PSICluster::getConnectivity() const {
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

void PSICluster::dumpCoefficients(std::ostream& os,
		PSICluster::ClusterPair const& curr) const {

	os << "a[0-2][0-2]: " << ' ' << curr.a00 << ' ' << curr.a01 << ' '
			<< curr.a02 << ' ' << curr.a10 << ' ' << curr.a11 << ' ' << curr.a12
			<< ' ' << curr.a20 << ' ' << curr.a21 << ' ' << curr.a22;
}

void PSICluster::dumpCoefficients(std::ostream& os,
		PSICluster::CombiningCluster const& curr) const {

	os << "a[0-2]: " << curr.a0 << ' ' << curr.a1 << ' ' << curr.a2;
}

void PSICluster::outputCoefficientsTo(std::ostream& os) const {

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

