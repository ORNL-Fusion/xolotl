#ifndef ALLOYCASES_H
#define ALLOYCASES_H

// Includes
#include "AlloyReactions.h"
#include <Constants.h>

namespace xolotlCore {

std::vector<forwardReaction> getForwardReactions (std::string caseName)
{
  std::vector<forwardReaction> reactions;
  forwardReaction reaction ("Cat","Dog");

  //if (caseName == "insitu") {
    // vac + vac = vac | void
    reaction = forwardReaction(vType,vType);
    reaction.addProduct(vType);
    reaction.addProduct(voidType);
    reactions.push_back(reaction);

    // vac + void = void
    reaction = forwardReaction(vType,voidType);
    reaction.addProduct(voidType);
    reactions.push_back(reaction);

    // vac + faulted = faulted
    reaction = forwardReaction(vType,faultedType);
    reaction.addProduct(faultedType);
    reactions.push_back(reaction);

    // vac + int = vac | int | recombine
    reaction = forwardReaction(vType,iType);
    reaction.addProduct(vType);
    reaction.addProduct(iType);
    reaction.addProduct("recombine");
    reactions.push_back(reaction);

    // vac + frank = vac | int | frank | recombine
    reaction = forwardReaction(vType,frankType);
    reaction.addProduct(vType);
    reaction.addProduct(iType);
    reaction.addProduct(frankType);
    reaction.addProduct("recombine");
    reactions.push_back(reaction);

    // vac + perfect = vac | int | perfect | recombine
    reaction = forwardReaction(vType,perfectType);
    reaction.addProduct(vType);
    reaction.addProduct(iType);
    reaction.addProduct(perfectType);
    reaction.addProduct("recombine");
    reactions.push_back(reaction);

    ///////////////////////////////////////////////////

    // int + void = void | vac | int | recombine
    reaction = forwardReaction(iType,voidType);
    reaction.addProduct(voidType);
    reaction.addProduct(vType);
    reaction.addProduct(iType);
    reaction.addProduct("recombine");
    reactions.push_back(reaction);

    // int + faulted = faulted | vac | int | recombine
    reaction = forwardReaction(iType,faultedType);
    reaction.addProduct(faultedType);
    reaction.addProduct(vType);
    reaction.addProduct(iType);
    reaction.addProduct("recombine");
    reactions.push_back(reaction);

    // int + int = int | frank
    reaction = forwardReaction(iType,iType);
    reaction.addProduct(iType);
    reaction.addProduct(frankType);
    reactions.push_back(reaction);

    // int + frank = frank
    reaction = forwardReaction(iType,frankType);
    reaction.addProduct(frankType);
    reactions.push_back(reaction);

    // int + perfect = perfect
    reaction = forwardReaction(iType,perfectType);
    reaction.addProduct(perfectType);
    reactions.push_back(reaction);

    ///////////////////////////////////////////////////

    // perfect + void = void | vac | int | perfect | recombine
    reaction = forwardReaction(perfectType,voidType);
    reaction.addProduct(voidType);
    reaction.addProduct(vType);
    reaction.addProduct(iType);
    reaction.addProduct(perfectType);
    reaction.addProduct("recombine");
    reactions.push_back(reaction);

    // perfect + faulted = faulted | vac | int | perfect | recombine
    reaction = forwardReaction(perfectType,faultedType);
    reaction.addProduct(faultedType);
    reaction.addProduct(vType);
    reaction.addProduct(iType);
    reaction.addProduct(perfectType);
    reaction.addProduct("recombine");
    reactions.push_back(reaction);

    // perfect + frank = frank
    reaction = forwardReaction(perfectType,frankType);
    reaction.addProduct(frankType);
    reactions.push_back(reaction);

    // perfect + perfect = perfect
    reaction = forwardReaction(perfectType,perfectType);
    reaction.addProduct(perfectType);
    reactions.push_back(reaction);

  //}

  return reactions;
}

std::vector<backwardReaction> getBackwardReactions (std::string caseName)
{
  std::vector<backwardReaction> reactions;
  backwardReaction reaction ("Cat","Dog");

  //if (caseName == "insitu") {
    // void = vac + (void | vac)
    reaction = backwardReaction(voidType,vType);
    reaction.addProduct(voidType);
    reaction.addProduct(vType);
    reactions.push_back(reaction);

    // vac = vac + (vac)
    reaction = backwardReaction(vType,vType);
    reaction.addProduct(vType);
    reactions.push_back(reaction);

    // faulted = vac + (faulted | vac)
    reaction = backwardReaction(faultedType,vType);
    reaction.addProduct(faultedType);
    reaction.addProduct(vType);
    reactions.push_back(reaction);


  //}

  return reactions;
}

} /* end namespace xolotlCore */
#endif
