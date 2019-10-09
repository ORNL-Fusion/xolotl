// Includes
#include <HeterogeneousNucleationHandler.h>
#include <Constants.h>

namespace xolotlCore {

void HeterogeneousNucleationHandler::initialize(
		const IReactionNetwork& network) {
	// Get the two smallest xenon clusters because they are the only ones involved
	auto singleXenon = network.get(Species::Xe, 1);
	auto doubleXenon = network.get(Species::Xe, 2);

	if (!singleXenon || !doubleXenon) {
		// Inform the user
		std::cout << "The heterogeneous nucleation won't happen because "
				"the single or double xenon cluster is missing." << std::endl;

		return;
	}

	// Add the connectivity
	singleXenon->setReactionConnectivity(doubleXenon->getId());
	doubleXenon->setReactionConnectivity(singleXenon->getId());

	return;
}

void HeterogeneousNucleationHandler::updateHeterogeneousNucleationRate(
		double rate) {
	// We say there are 25 bubbles created per fission fragments and there are 2 fission fragments per fission
	nucleationRate = 50.0 * rate;

	return;
}

void HeterogeneousNucleationHandler::computeHeterogeneousNucleation(
		const IReactionNetwork& network, double *concOffset,
		double *updatedConcOffset, int xi, int xs, int yj, int zk) {
	// Get the single and double xenon
	auto singleXenon = network.get(Species::Xe, 1), doubleXenon = network.get(
			Species::Xe, 2);
	int singleXenonId = singleXenon->getId() - 1, doubleXenonId =
			doubleXenon->getId() - 1;

	// Get the single concentration to know in which regime we are
	double singleConc = singleXenon->getConcentration();
	// Update the concentrations
	if (singleConc > 2.0 * nucleationRate) {
		updatedConcOffset[singleXenonId] -= 2.0 * nucleationRate;
		updatedConcOffset[doubleXenonId] += nucleationRate;
	} else {
		updatedConcOffset[singleXenonId] -= singleConc;
		updatedConcOffset[doubleXenonId] += singleConc / 2.0;
	}

	// Remove the contribution from homogeneous nucleation
	updatedConcOffset[singleXenonId] += 32.0 * xolotlCore::pi
			* singleXenon->getReactionRadius()
			* singleXenon->getDiffusionCoefficient(0) * singleConc * singleConc;
	updatedConcOffset[doubleXenonId] -= 16.0 * xolotlCore::pi
			* singleXenon->getReactionRadius()
			* singleXenon->getDiffusionCoefficient(0) * singleConc * singleConc;

	return;
}

void HeterogeneousNucleationHandler::computePartialsForHeterogeneousNucleation(
		const IReactionNetwork& network, double *val, int *indices, int xi,
		int xs, int yj, int zk) {
	// Get the single and double xenon
	auto singleXenon = network.get(Species::Xe, 1), doubleXenon = network.get(
			Species::Xe, 2);
	int singleXenonId = singleXenon->getId() - 1, doubleXenonId =
			doubleXenon->getId() - 1;

	// Get the single concentration to know in which regime we are
	double singleConc = singleXenon->getConcentration();
	// Set the indices
	indices[0] = singleXenonId;
	indices[1] = doubleXenonId;
	// Update the partials
	if (singleConc > 2.0 * nucleationRate) {
		val[0] = 0.0;
		val[1] = 0.0;
	} else {
		val[0] = -1.0;
		val[1] = 0.5;
	}
	// Remove the contribution from homogeneous nucleation
	val[0] += 32.0 * xolotlCore::pi * singleXenon->getReactionRadius()
			* singleXenon->getDiffusionCoefficient(0) * singleConc;
	val[1] -= 16.0 * xolotlCore::pi * singleXenon->getReactionRadius()
			* singleXenon->getDiffusionCoefficient(0) * singleConc;

	return;
}

}/* end namespace xolotlCore */

