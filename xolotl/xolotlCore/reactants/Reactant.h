#ifndef REACTANT_H
#define REACTANT_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include <map>
//#include "ReactionNetwork.h"

namespace xolotlCore {

class ReactionNetwork;

/**
 * A Reactant is a reacting body in a reaction network. It represents any body
 * whose population can change with time due to reactions of any type.
 */
class Reactant
{

    protected:

        //! The total concentration of this Reactant.
        double concentration;

        //! The name of this Reactant.
        std::string name;

        //! An integer identification number for this reactant.
        int id;

    public:

        //! The constructor.
        Reactant();

        //! Destructor
        ~Reactant(); 

        /**
         * An alternative constructor that can be used to create a reactant
         * with an initial concentration.
         * @param conc The initial concentration
         */
        Reactant(double conc); 

        /**
         * This operation returns the current concentration.
         * @return The concentration of this reactant
         */
        double getConcentration(); 

        /**
         * This operation increases the concentration of the reactant by the
         * specified amount.
         * @param deltaConc the change in concentration
         */
        void increaseConcentration(double deltaConc); 

        /**
         * This operation decreases the concentration of the reactant by the
         * specified amount.
         * @param deltaConc the change in concentration
         */
        void decreaseConcentration(double deltaConc); 

        /**
         * This operation sets the concentration of the reactant to the
         * specified amount.
         * @param conc The new concentation
         */
        void setConcentration(double conc); 

        /**
         * This operation sets the concentration of the reactant to zero.
         */
        void zero(); 

        /**
         * This operation returns the total flux of this reactant in the
         * current network.
         * @return The total change in flux for this reactant due to all
         * reactions
         */
        virtual double getTotalFlux();

        /**
         * This operation sets the collection of other reactants that make up
         * the reaction network in which this reactant exists.
         * @param network The reaction network of which this reactant is a part
         */
        void setReactionNetwork(std::shared_ptr<ReactionNetwork> network);

        /**
         * This operation writes the contents of the reactant to a string. This
         * operation should be overridden by subclasses.
         * @return A serialized version of this reactant as a string.
         */
        virtual const std::string toString();

};  //end class Reactant

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

} /* namespace xolotlCore */
#endif
