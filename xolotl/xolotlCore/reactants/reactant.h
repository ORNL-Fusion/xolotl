#ifndef REACTANT_H
#define REACTANT_H

// Includes
#include <string>
#include <vector>

namespace xolotlCore {

/**
 * A Reactant is a reacting body in a reaction network. It represents any body
 * whose population can change with time due to reactions of any type.
 */
class Reactant
{

    private:

        // The total concentration of this Reactant.
        double concentration;

        // The name of this Reactant. 
        const std::string name;

        // An integer identification number for this reactant.  
        int id;

    public:

        // The constructor.
        Reactant();

        // Destructor
        ~Reactant(); 

        // An alternative constructor that can be used to create a reactant with an initial concentration. 
        Reactant(double conc); 

        // This operation returns the current concentration. 
        double getConcentration(); 

        // This operation increases the concentration of the reactant by the specified amount. 
        void increaseConcentration(double deltaConc); 

        // This operation decreases the concentration of the reactant by the specified amount. 
        void decreaseConcentration(double deltaConc); 

        // This operation sets the concentration of the reactant to the specified amount. 
        void setConcentration(double conc); 

        // This operation sets the concentration of the reactant to zero. 
        void zero(); 

        // This operation returns the total flux of this reactant in the current network. 
        double getTotalFlux(); 

        // This operation sets the collection of other reactants that make up the reaction network in which this reactant exists.  
        void setReactionNetwork(std::vector<Reactant> network);

        // This operation writes the contents of the reactant to a string. 
        const std::string toString();

        // This operation reads the contents of the reactant from a string. 
        void fromString(const std::string input);

        // This operation writes the contents of the reactant to a string in XML. 
        const std::string toXML();

        // This operation reads the contents of the reactant from a string in XML. 
        void fromXML(const std::string input);

};  //end class Reactant

} /* namespace xolotlCore */
#endif
