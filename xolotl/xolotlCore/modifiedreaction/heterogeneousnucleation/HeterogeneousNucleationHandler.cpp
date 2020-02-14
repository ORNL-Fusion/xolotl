// Includes
#include <HeterogeneousNucleationHandler.h>
#include <experimental/NEReactionNetwork.h>
#include <Constants.h>

namespace xolotlCore {

void HeterogeneousNucleationHandler::initialize(
		experimental::IReactionNetwork& network,
		experimental::IReactionNetwork::SparseFillMap& dfill) {

	using NetworkType =
	experimental::NEReactionNetwork;
	auto neNetwork = dynamic_cast<NetworkType*>(&network);
	// Get the two smallest xenon clusters because they are the only ones involved
	NetworkType::Composition comp;
	for (auto i : neNetwork->getSpeciesRange()) {
		comp[i] = 0;
	}
	comp[NetworkType::Species::Xe] = 1;
	auto singleXenon = neNetwork->findCluster(comp, plsm::onHost);
	auto singleId = singleXenon.getId();
	comp[NetworkType::Species::Xe] = 2;
	auto doubleXenon = neNetwork->findCluster(comp, plsm::onHost);
	auto doubleId = doubleXenon.getId();

	if (singleId == plsm::invalid<std::size_t>
			|| doubleId == plsm::invalid<std::size_t>) {
		// Inform the user
		std::cout << "The heterogeneous nucleation won't happen because "
				"the single or double xenon cluster is missing." << std::endl;

		return;
	}

	// Add the connectivity
	dfill[singleId].emplace_back(singleId);
	dfill[singleId].emplace_back(doubleId);
	dfill[doubleId].emplace_back(singleId);
	dfill[doubleId].emplace_back(doubleId);

	return;
}

void HeterogeneousNucleationHandler::updateHeterogeneousNucleationRate(
		double rate) {
	// We say there are 25 bubbles created per fission fragments and there are 2 fission fragments per fission
	nucleationRate = 50.0 * rate / fissionYield;

	return;
}

void HeterogeneousNucleationHandler::setFissionYield(double yield) {
	fissionYield = yield;

	return;
}

void HeterogeneousNucleationHandler::computeHeterogeneousNucleation(
		experimental::IReactionNetwork& network, double *concOffset,
		double *updatedConcOffset, int xi, int xs, int yj, int zk) {

	// TODO: it might be interesting to save this information instead of getting it every time

	using NetworkType =
	experimental::NEReactionNetwork;
	auto neNetwork = dynamic_cast<NetworkType*>(&network);
	// Get the two smallest xenon clusters because they are the only ones involved
	NetworkType::Composition comp;
	for (auto i : neNetwork->getSpeciesRange()) {
		comp[i] = 0;
	}
	comp[NetworkType::Species::Xe] = 1;
	auto singleXenon = neNetwork->findCluster(comp, plsm::onHost);
	auto singleId = singleXenon.getId();
	comp[NetworkType::Species::Xe] = 2;
	auto doubleXenon = neNetwork->findCluster(comp, plsm::onHost);
	auto doubleId = doubleXenon.getId();

	// Get the single concentration to know in which regime we are
	double singleConc = concOffset[singleId];

	// Update the concentrations
	if (singleConc > 2.0 * nucleationRate) {
		updatedConcOffset[singleId] -= 2.0 * nucleationRate;
		updatedConcOffset[doubleId] += nucleationRate;
	} else {
		updatedConcOffset[singleId] -= singleConc;
		updatedConcOffset[doubleId] += singleConc / 2.0;
	}

	// Remove the contribution from homogeneous nucleation
	updatedConcOffset[singleId] += 32.0 * xolotlCore::pi
			* singleXenon.getReactionRadius()
			* singleXenon.getDiffusionCoefficient(0) * singleConc * singleConc;
	updatedConcOffset[doubleId] -= 16.0 * xolotlCore::pi
			* singleXenon.getReactionRadius()
			* singleXenon.getDiffusionCoefficient(0) * singleConc * singleConc;

	return;
}

bool HeterogeneousNucleationHandler::computePartialsForHeterogeneousNucleation(
		experimental::IReactionNetwork& network, double *concOffset,
		double *val, int *indices, int xi, int xs, int yj, int zk) {

	using NetworkType =
	experimental::NEReactionNetwork;
	auto neNetwork = dynamic_cast<NetworkType*>(&network);
	// Get the two smallest xenon clusters because they are the only ones involved
	NetworkType::Composition comp;
	for (auto i : neNetwork->getSpeciesRange()) {
		comp[i] = 0;
	}
	comp[NetworkType::Species::Xe] = 1;
	auto singleXenon = neNetwork->findCluster(comp, plsm::onHost);
	auto singleId = singleXenon.getId();
	comp[NetworkType::Species::Xe] = 2;
	auto doubleXenon = neNetwork->findCluster(comp, plsm::onHost);
	auto doubleId = doubleXenon.getId();

	// Get the single concentration to know in which regime we are
	double singleConc = concOffset[singleId];

	// Set the indices
	indices[0] = singleId;
	indices[1] = doubleId;
	// Update the partials
	if (singleConc > 2.0 * nucleationRate) {
		val[0] = 0.0;
		val[1] = 0.0;
	} else {
		val[0] = -1.0;
		val[1] = 0.5;
	}
	// Remove the contribution from homogeneous nucleation
	val[0] += 32.0 * xolotlCore::pi * singleXenon.getReactionRadius()
			* singleXenon.getDiffusionCoefficient(0) * singleConc;
	val[1] -= 16.0 * xolotlCore::pi * singleXenon.getReactionRadius()
			* singleXenon.getDiffusionCoefficient(0) * singleConc;

	return true;
}

}/* end namespace xolotlCore */

