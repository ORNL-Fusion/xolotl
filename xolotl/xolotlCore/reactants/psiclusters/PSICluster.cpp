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
	double firstDistance[3] = { 1.0, 0.0, 0.0 }, secondDistance[3] = { 1.0, 0.0,
			0.0 };
	if (newPair.first.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(newPair.first);
		firstDistance[1] = super.getDistance(c, 0);
		firstDistance[2] = super.getDistance(d, 1);
	}
	if (newPair.second.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(newPair.second);
		secondDistance[1] = super.getDistance(c, 0);
		secondDistance[2] = super.getDistance(d, 1);
	}
	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < 3; i++) {
			newPair.coefs[i][j] += firstDistance[i] * secondDistance[j];
		}
	}

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

				double firstDistance[3] = {1.0, 0.0, 0.0}, secondDistance[3] = {1.0, 0.0, 0.0};
				if (newPair.first.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster const&>(newPair.first);
					firstDistance[1] = super.getDistance(c, 0);
					firstDistance[2] = super.getDistance(d, 1);
				}
				if (newPair.second.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSICluster const&>(newPair.second);
					secondDistance[1] = super.getDistance(c, 0);
					secondDistance[2] = super.getDistance(d, 1);
				}
				for (int j = 0; j < 3; j++) {
					for (int i = 0; i < 3; i++) {
						newPair.coefs[i][j] += firstDistance[i] * secondDistance[j];
					}
				}
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
		heDistance = otherCluster.getDistance(a, 0);
		vDistance = otherCluster.getDistance(b, 1);
	}
	(*it).coefs[0] += 1.0;
	(*it).coefs[1] += heDistance;
	(*it).coefs[2] += vDistance;

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
					heDistance = otherCluster.getDistance(a, 0);
					vDistance = otherCluster.getDistance(b, 1);
				}
				combCluster.coefs[0] += 1.0;
				combCluster.coefs[1] += heDistance;
				combCluster.coefs[2] += vDistance;
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
		firstHeDistance = super.getDistance(a, 0);
		firstVDistance = super.getDistance(b, 1);
	}
	(*it).coefs[0][0] += 1.0;
	(*it).coefs[1][0] += firstHeDistance;
	(*it).coefs[2][0] += firstVDistance;

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
					firstHeDistance = super.getDistance(a, 0);
					firstVDistance = super.getDistance(b, 1);
				}
				currPair.coefs[0][0] += 1.0;
				currPair.coefs[1][0] += firstHeDistance;
				currPair.coefs[2][0] += firstVDistance;
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
				double lA[3] = {};
				lA[0] = dissCluster.getConcentration(0.0, 0.0);
				lA[1] = dissCluster.getHeMomentum();
				lA[2] = dissCluster.getVMomentum();

				double sum = 0.0;
				for (int i = 0; i < 3; i++) {
					sum += currPair.coefs[i][0] * lA[i];
				}

				// Calculate the Dissociation flux
				return running +
				(currPair.reaction.kConstant * sum);
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
			double lA[3] = {}, lB[3] = {};
			lA[0] = firstReactant.getConcentration(0.0, 0.0);
			lB[0] = secondReactant.getConcentration(0.0, 0.0);
			lA[1] = firstReactant.getHeMomentum();
			lB[1] = secondReactant.getHeMomentum();
			lA[2] = firstReactant.getVMomentum();
			lB[2] = secondReactant.getVMomentum();

			double sum = 0.0;
			for (int j = 0; j < 3; j++) {
				for (int i = 0; i < 3; i++) {
					sum += currPair.coefs[i][j] * lA[i] * lB[j];
				}
			}
			// Update the flux
			return running + (currPair.reaction.kConstant *
					sum);
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
				double lB[3] = {};
				lB[0] = combiningCluster.getConcentration(0.0, 0.0);
				lB[1] = combiningCluster.getHeMomentum();
				lB[2] = combiningCluster.getVMomentum();

				double sum = 0.0;
				for (int i = 0; i < 3; i++) {
					sum += cc.coefs[i] * lB[i];
				}
				// Calculate the combination flux
				return running + (cc.reaction.kConstant *
						sum);

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
				double lA[3] = {}, lB[3] = {};
				lA[0] = firstReactant.getConcentration(0.0, 0.0);
				lB[0] = secondReactant.getConcentration(0.0, 0.0);
				lA[1] = firstReactant.getHeMomentum();
				lB[1] = secondReactant.getHeMomentum();
				lA[2] = firstReactant.getVMomentum();
				lB[2] = secondReactant.getVMomentum();

				// Compute contribution from the first part of the reacting pair
				double value = currPair.reaction.kConstant;

				double sum[3][2] = {};
				for (int j = 0; j < 3; j++) {
					for (int i = 0; i < 3; i++) {
						sum[j][0] += currPair.coefs[j][i] * lB[i];
						sum[j][1] += currPair.coefs[i][j] * lA[i];
					}
				}

				partials[firstReactant.id - 1] += value * sum[0][0];
				partials[firstReactant.heMomId - 1] += value * sum[1][0];
				partials[firstReactant.vMomId - 1] += value * sum[2][0];

				// Compute contribution from the second part of the reacting pair
				partials[secondReactant.id - 1] += value * sum[0][1];
				partials[secondReactant.heMomId - 1] += value * sum[1][1];
				partials[secondReactant.vMomId - 1] += value * sum[2][1];
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
				* (cc.coefs[0] * l0B + cc.coefs[1] * lHeB + cc.coefs[2] * lVB);
				// Compute the contribution from the combining cluster
				double value = cc.reaction.kConstant * concentration;
				partials[cluster.id - 1] -= value * cc.coefs[0];
				partials[cluster.heMomId - 1] -= value * cc.coefs[1];
				partials[cluster.vMomId - 1] -= value * cc.coefs[2];
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
				partials[cluster.id - 1] += value * currPair.coefs[0][0];
				partials[cluster.heMomId - 1] += value * currPair.coefs[1][0];
				partials[cluster.vMomId - 1] += value * currPair.coefs[2][0];
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
						(cc.reaction.kConstant * cc.combining.concentration * cc.coefs[0]);
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

	os << "a[0-2][0-2]: " << ' ' << curr.coefs[0][0] << ' ' << curr.coefs[0][1]
			<< ' ' << curr.coefs[0][2] << ' ' << curr.coefs[1][0] << ' '
			<< curr.coefs[1][1] << ' ' << curr.coefs[1][2] << ' '
			<< curr.coefs[2][0] << ' ' << curr.coefs[2][1] << ' '
			<< curr.coefs[2][2];
}

void PSICluster::dumpCoefficients(std::ostream& os,
		PSICluster::CombiningCluster const& curr) const {

	os << "a[0-2]: " << curr.coefs[0] << ' ' << curr.coefs[1] << ' '
			<< curr.coefs[2];
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

