// Includes
#include "Reactant.h"
#include "ReactionNetwork.h"
#include <xolotlPerf.h>
#include <iostream>
#include <math.h>
#include <MathUtils.h>
#include <Constants.h>

namespace xolotlCore {

// TODO modify signature to take type as argument.
Reactant::Reactant(IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry,
		const std::string& _name) :
		concentration(0.0), id(0), temperature(0.0), type(
				ReactantType::Invalid), network(_network), handlerRegistry(
				registry), size(0), formationEnergy(0.0), diffusionFactor(0.0), diffusionCoefficient(
				0.0), migrationEnergy(0.0), name(_name), reactionRadius(0.0) {

}

void Reactant::recomputeDiffusionCoefficient(double temp) {
	// Return zero if the diffusion factor is zero.
	if (xolotlCore::equal(diffusionFactor, 0.0)) {
		diffusionCoefficient = 0.0;
	} else {
		// Otherwise use the Arrhenius equation to compute the diffusion
		// coefficient
		double k_b = xolotlCore::kBoltzmann;
		double kernel = -1.0 * migrationEnergy / (k_b * temp);
		diffusionCoefficient = diffusionFactor * exp(kernel);
	}

	return;
}

std::vector<int> Reactant::getConnectivity() const {
	// The connectivity array by default is filled with
	// zeros.
	int connectivityLength = network.getDOF();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);

	// This reactant should be connected to itself
	connectivity[id - 1] = 1;

	return connectivity;
}

void Reactant::setTemperature(double temp) {
	temperature = temp;

	// Recompute the diffusion coefficient
	recomputeDiffusionCoefficient(temp);
}

void Reactant::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

void Reactant::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

std::ostream&
operator<<(std::ostream& os, const IReactant::Composition& comp) {
	std::vector<Species> compSpecies { Species::He, Species::I, Species::V,
			Species::Xe };
	for (auto const& currSpecies : compSpecies) {
		os << toString(currSpecies) << comp[toCompIdx(currSpecies)];
	}
	return os;
}

std::ostream&
operator<<(std::ostream& os, const IReactant& reactant) {
	os << "id: " << reactant.getId() << "; " << "type: "
			<< toString(reactant.getType()) << "; " << "comp: "
			<< reactant.getComposition();
	return os;
}

} // namespace xolotlCore
