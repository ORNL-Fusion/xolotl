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
		concentration(0.0), name(_name), type(ReactantType::Invalid), id(0), network(
				_network), handlerRegistry(registry), size(0), diffusionFactor(
				0.0), formationEnergy(0.0), migrationEnergy(
				0.0), reactionRadius(0.0) {
}

void Reactant::recomputeDiffusionCoefficient(double temp, int i) {
	// Return zero if the diffusion factor is zero.
	if (!xolotlCore::equal(diffusionFactor, 0.0)) {
		// Otherwise use the Arrhenius equation to compute the diffusion
		// coefficient
		double k_b = xolotlCore::kBoltzmann;
		double kernel = -1.0 * migrationEnergy / (k_b * temp);
		diffusionCoefficient[i] = diffusionFactor * exp(kernel);
	}

	return;
}

void Reactant::addGridPoints(int i) {
	// Add grid points
	if (i > 0) {
		while (i > 0) {
			diffusionCoefficient.emplace(diffusionCoefficient.begin(), 0.0);
			temperature.emplace(temperature.begin(), 0.0);

			// Decrease i
			i--;
		}
	} else {
		diffusionCoefficient.erase(diffusionCoefficient.begin(),
				diffusionCoefficient.begin() - i);
		temperature.erase(temperature.begin(), temperature.begin() - i);
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

void Reactant::setTemperature(double temp, int i) {
	temperature[i] = temp;

	// Recompute the diffusion coefficient
	recomputeDiffusionCoefficient(temp, i);
}

void Reactant::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;

	return;
}

void Reactant::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;

	return;
}

std::ostream&
operator<<(std::ostream& os, const IReactant::Composition& comp) {
	std::vector<Species> compSpecies { Species::He, Species::D, Species::T,
			Species::I, Species::V, Species::Xe };
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
