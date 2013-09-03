#include "ReactionNetwork.h"
#include "Reactant.h"

using namespace xolotlCore;

ReactionNetwork::ReactionNetwork(const ReactionNetwork &other) {
	// The copy constructor of std::map copies each of the keys and values.
	properties.reset(new std::map<std::string, std::string>(*other.properties));
}

std::map<std::string, int> ReactionNetwork::toClusterMap(int index) const {
	// This base class returns an empty map
	return std::map<std::string, int>();
}

int ReactionNetwork::toClusterIndex(std::map<std::string, int> clusterMap) const {
	// This base class returns a zero index value
	return 0;
}

std::map<std::string,int> Reactant::getComposition() const {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Return an empty map
	return clusterMap;
}
