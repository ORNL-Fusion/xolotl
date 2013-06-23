
// Includes
#include "Reactant.h"

// Namespaces
using namespace xolotlCore;

Reactant::Reactant() :
concentration(0.0),name(""),id(0)
{
    //TODO Auto-generated method stub
}
Reactant::~Reactant() 
{
    //TODO Auto-generated method stub
}
Reactant::Reactant(double conc) : 
concentration(0.0),name(""),id(0)
{
    //TODO Auto-generated method stub
}
double Reactant::getConcentration() 
{
    //TODO Auto-generated method stub
    return 0;
}
void Reactant::increaseConcentration(double deltaConc) 
{
    //TODO Auto-generated method stub
    return;
}
void Reactant::decreaseConcentration(double deltaConc) 
{
    //TODO Auto-generated method stub
    return;
}
void Reactant::setConcentration(double conc) 
{
    //TODO Auto-generated method stub
    return;
}
void Reactant::zero() 
{
    //TODO Auto-generated method stub
    return;
}
double Reactant::getTotalFlux() 
{
    //TODO Auto-generated method stub
    return 0;
}
void Reactant::setReactionNetwork(std::shared_ptr<ReactionNetwork> network)
{
    //TODO Auto-generated method stub
    return;
}

std::vector<int> Reactant::getConnectivity() {
	std::vector<int> dummy;
	return dummy;
}

const std::string Reactant::toString()
{
    //TODO Auto-generated method stub
    return "fluffy";
}

