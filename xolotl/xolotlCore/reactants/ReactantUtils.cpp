#include <sstream>
#include "ReactionNetwork.h"

namespace xolotlCore {

ProductionReaction::ProductionReaction(IReactant * _reactant1, IReactant * _reactant2)
  : first(_reactant1),
    second(_reactant2) {

    // We made an assumption about ordering.
    // Check if we were wrong, and if so, swap.
    if(_reactant2->getCompositionString() < _reactant1->getCompositionString()) {
        auto tmp = first;
        first = second;
        second = tmp;
    }

    // Determine our descriptive key and cache it.
    // This assumes that our reactants are ordered.
    descKey = first->getCompositionString() + second->getCompositionString();
}

DissociationReaction::DissociationReaction(IReactant * dissociatingPtr, IReactant * _reactant1, IReactant * _reactant2)
  : dissociating(dissociatingPtr),
    first(_reactant1),
    second(_reactant2),
    reverseReaction(nullptr) {

    // We made an assumption about ordering.
    // Check if we were wrong, and if so, swap.
    if(_reactant2->getCompositionString() < _reactant1->getCompositionString()) {
        auto tmp = first;
        first = second;
        second = tmp;
    }

    // Determine our descriptive key and cache it.
    // This assumes that our reactants are ordered.
    descKey = dissociating->getCompositionString()
        + first->getCompositionString() 
        + second->getCompositionString();
}


} // namespace xolotlCore

