#include "ReactionNetwork.h"
#include "Reactant.h"

using namespace xolotlCore;

ReactionNetwork::ReactionNetwork(const ReactionNetwork &other) {
	// The copy constructor of std::map copies each of the keys and values.
	properties.reset(new std::map<std::string, std::string>(*other.properties));

	// Copy the reactants list by pushing a copy of each element onto the
	// reactants vector
	reactants.reset(new std::vector<std::shared_ptr<Reactant>>);

	for (std::vector<std::shared_ptr<Reactant>>::iterator reactantIt =
			other.reactants->begin(); reactantIt != other.reactants->end();
			reactantIt++) {

		std::shared_ptr<Reactant> reactant(new Reactant(**reactantIt));
		reactants->push_back(reactant);
	}
}

std::map<std::string, int> ReactionNetwork::toClusterMap(int index) const {
	// This base class returns an empty map
	return std::map<std::string, int>();
}

int ReactionNetwork::toClusterIndex(std::map<std::string, int> clusterMap) const {
	// This base class returns a zero index value
	return 0;
}
