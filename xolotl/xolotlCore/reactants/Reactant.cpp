// Includes
#include "Reactant.h"
#include "ReactionNetwork.h"
#include <xolotlPerf.h>
#include <iostream>
#include <math.h>
#include <MathUtils.h>
#include <Constants.h>

// Namespaces
using namespace xolotlCore;

Reactant::Reactant() :
		concentration(0.0), id(0), momId(0), heMomId(0), vMomId(0), temperature(
				0.0), typeName(""), network(nullptr), handlerRegistry(nullptr), size(
				0), formationEnergy(0.0), diffusionFactor(0.0), diffusionCoefficient(
				0.0), migrationEnergy(0.0), name("Reactant"), reactionRadius(
				0.0) {
}

Reactant::Reactant(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		concentration(0.0), id(0), momId(0), heMomId(0), vMomId(0), temperature(
				0.0), typeName(""), network(nullptr), handlerRegistry(registry), size(
				0), formationEnergy(0.0), diffusionFactor(0.0), diffusionCoefficient(
				0.0), migrationEnergy(0.0), name("Reactant"), reactionRadius(
				0.0) {
}

Reactant::Reactant(Reactant &other) :
		concentration(other.concentration), name(other.name), typeName(
				other.typeName), id(other.id), momId(other.momId), heMomId(
				other.heMomId), vMomId(other.vMomId), temperature(
				other.temperature), network(other.network), handlerRegistry(
				other.handlerRegistry), size(other.size), formationEnergy(
				other.formationEnergy), diffusionFactor(other.diffusionFactor), diffusionCoefficient(
				other.diffusionCoefficient), migrationEnergy(
				other.migrationEnergy), reactionRadius(other.reactionRadius), reactionConnectivitySet(
				other.reactionConnectivitySet), dissociationConnectivitySet(
				other.dissociationConnectivitySet) {
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
	int connectivityLength = network->getDOF();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);

	// This reactant should be connected to itself
	connectivity[id - 1] = 1;

	return connectivity;
}

std::string Reactant::toCanonicalString(std::string type,
		const std::map<std::string, int>& composition) {

	// Construct the canonical string representation of the given composition.
	// Note that this can only be considered to produce a canonical
	// representation of the reactant map's contents if the reactant map
	// is kept in a collection where we are guaranteed that iterating
	// over its contents will always be done in the same order, regardless
	// of how the map was produced.  Thankfully, std::map is an
	// ordered map and gives this guarantee.
	std::ostringstream ostr;
	ostr << type << ':';
	for (auto iter = composition.begin(); iter != composition.end(); ++iter) {
		// Add the current reactant's name and size to the string.
		// Note that we don't really care about nice formatting, since
		// this isn't intended to be a human-readable string.
		ostr << iter->first << iter->second;
	}
	return ostr.str();
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

