#ifndef ALLOYCASES_H
#define ALLOYCASES_H

// Includes
#include "AlloyReactions.h"
#include <Constants.h>

namespace xolotlCore {

std::vector<forwardReaction> getForwardReactions(std::string caseName) {
	std::vector<forwardReaction> reactions;

	//if (caseName == "insitu") {
	// vac + vac = vac | void
	forwardReaction reaction = forwardReaction(ReactantType::V,
			ReactantType::V);
	reaction.addProduct(ReactantType::V);
	reaction.addProduct(ReactantType::Void);
	reactions.push_back(reaction);

	// vac + void = void
	reaction = forwardReaction(ReactantType::V, ReactantType::Void);
	reaction.addProduct(ReactantType::Void);
	reactions.push_back(reaction);

	// vac + faulted = faulted
	reaction = forwardReaction(ReactantType::V, ReactantType::Faulted);
	reaction.addProduct(ReactantType::Faulted);
	reactions.push_back(reaction);

	// vac + int = vac | int | recombine
	reaction = forwardReaction(ReactantType::V, ReactantType::I);
	reaction.addProduct(ReactantType::V);
	reaction.addProduct(ReactantType::I);
	reactions.push_back(reaction);

	// vac + frank = vac | int | frank | recombine
	reaction = forwardReaction(ReactantType::V, ReactantType::Frank);
	reaction.addProduct(ReactantType::V);
	reaction.addProduct(ReactantType::I);
	reaction.addProduct(ReactantType::Frank);
	reactions.push_back(reaction);

	// vac + perfect = vac | int | perfect | recombine
	reaction = forwardReaction(ReactantType::V, ReactantType::Perfect);
	reaction.addProduct(ReactantType::V);
	reaction.addProduct(ReactantType::I);
	reaction.addProduct(ReactantType::Perfect);
	reactions.push_back(reaction);

	///////////////////////////////////////////////////

	// int + void = void | vac | int | recombine
	reaction = forwardReaction(ReactantType::I, ReactantType::Void);
	reaction.addProduct(ReactantType::Void);
	reaction.addProduct(ReactantType::V);
	reaction.addProduct(ReactantType::I);
	reactions.push_back(reaction);

	// int + faulted = faulted | vac | int | recombine
	reaction = forwardReaction(ReactantType::I, ReactantType::Faulted);
	reaction.addProduct(ReactantType::Faulted);
	reaction.addProduct(ReactantType::V);
	reaction.addProduct(ReactantType::I);
	reactions.push_back(reaction);

	// int + int = int | frank
	reaction = forwardReaction(ReactantType::I, ReactantType::I);
	reaction.addProduct(ReactantType::I);
	reaction.addProduct(ReactantType::Frank);
	reactions.push_back(reaction);

	// int + frank = frank
	reaction = forwardReaction(ReactantType::I, ReactantType::Frank);
	reaction.addProduct(ReactantType::Frank);
	reactions.push_back(reaction);

	// int + perfect = perfect
	reaction = forwardReaction(ReactantType::I, ReactantType::Perfect);
	reaction.addProduct(ReactantType::Perfect);
	reactions.push_back(reaction);

	///////////////////////////////////////////////////

	// perfect + void = void | vac | int | perfect | recombine
	reaction = forwardReaction(ReactantType::Perfect, ReactantType::Void);
	reaction.addProduct(ReactantType::Void);
	reaction.addProduct(ReactantType::V);
	reaction.addProduct(ReactantType::I);
	reaction.addProduct(ReactantType::Perfect);
	reactions.push_back(reaction);

	// perfect + faulted = faulted | vac | int | perfect | recombine
	reaction = forwardReaction(ReactantType::Perfect, ReactantType::Faulted);
	reaction.addProduct(ReactantType::Faulted);
	reaction.addProduct(ReactantType::V);
	reaction.addProduct(ReactantType::I);
	reaction.addProduct(ReactantType::Perfect);
	reactions.push_back(reaction);

	// perfect + frank = frank
	reaction = forwardReaction(ReactantType::Perfect, ReactantType::Frank);
	reaction.addProduct(ReactantType::Frank);
	reactions.push_back(reaction);

	// perfect + perfect = perfect
	reaction = forwardReaction(ReactantType::Perfect, ReactantType::Perfect);
	reaction.addProduct(ReactantType::Perfect);
	reactions.push_back(reaction);

	//}

	return reactions;
}

std::vector<backwardReaction> getBackwardReactions(std::string caseName) {
	std::vector<backwardReaction> reactions;

//	// void = vac + (void | vac)
//	backwardReaction reaction = backwardReaction(ReactantType::Void, ReactantType::V);
//	reaction.addProduct(ReactantType::Void);
//	reaction.addProduct(ReactantType::V);
//	reactions.push_back(reaction);
//
//	// vac = vac + (vac)
//	reaction = backwardReaction(ReactantType::V, ReactantType::V);
//	reaction.addProduct(ReactantType::V);
//	reactions.push_back(reaction);
//
//	// faulted = vac + (faulted | vac)
//	reaction = backwardReaction(ReactantType::Faulted, ReactantType::V);
//	reaction.addProduct(ReactantType::Faulted);
//	reaction.addProduct(ReactantType::V);
//	reactions.push_back(reaction);

	// bubble -> int
	backwardReaction reaction = backwardReaction(ReactantType::Void,
			ReactantType::I);
	reaction.addProduct(ReactantType::Void);
	reactions.push_back(reaction);

	// bubble -> vac
	reaction = backwardReaction(ReactantType::Void, ReactantType::V);
	reaction.addProduct(ReactantType::Void);
	reaction.addProduct(ReactantType::V);
	reactions.push_back(reaction);

	// faulted -> vac
	reaction = backwardReaction(ReactantType::Faulted, ReactantType::V);
	reaction.addProduct(ReactantType::Faulted);
	reaction.addProduct(ReactantType::V);
	reactions.push_back(reaction);

	// int -> int
	reaction = backwardReaction(ReactantType::I, ReactantType::I);
	reaction.addProduct(ReactantType::I);
	reactions.push_back(reaction);

	return reactions;
}

} /* end namespace xolotlCore */
#endif
