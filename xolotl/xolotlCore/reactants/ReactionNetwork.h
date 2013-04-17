/*
 * ReactionNetwork.h
 *
 *  Created on: Apr 17, 2013
 *      Author: bkj
 */

#ifndef REACTIONNETWORK_H_
#define REACTIONNETWORK_H_

// Includes
#include "Reactant.h"
#include <map>
#include <string>
#include <vector>
#include <memory>

namespace xolotlCore {

/**
 *  This is a simple convenience class that contains a set of Reactants and the
 *  properties that describe that set.
 */
class ReactionNetwork {
public:

	/**
	 * The properties of this network. The exact configuration of the map is
	 * specified by the class that loaded the network.
	 */
	std::shared_ptr<std::map<std::string, std::string>> properties;

	/**
	 * The set of reactants. The exact order of the reactants in the set is
	 * specified by the class that loaded them.
	 */
	std::shared_ptr<std::vector<Reactant>> reactants;

	//! The constructor. It initializes the properties map and reactants vector.
	ReactionNetwork() :
			properties(new std::map<std::string, std::string>()), reactants(
					new std::vector<Reactant>()) {
	}

	//! The destructor
	virtual ~ReactionNetwork() {
	}
};

} /* end namespace xolotlCore */
#endif /* REACTIONNETWORK_H_ */
