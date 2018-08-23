#include "AlloyReactions.h"

using namespace xolotlCore;

forwardReaction::forwardReaction (std::string name1, std::string name2)
{
  reactant1 = name1;
  reactant2 = name2;
  products.resize(0);
  return;
}

void forwardReaction::addProduct (std::string name)
{
  products.push_back(name);
  return;
}

std::string forwardReaction::getFirstReactant () const
{
  return reactant1;
}

std::string forwardReaction::getSecondReactant () const
{
  return reactant2;
}

std::vector<std::string> forwardReaction::getProducts () const
{
  return products;
}

backwardReaction::backwardReaction (std::string parentName,
    std::string monomerName)
{
  parent = parentName;
  monomer = monomerName;
  products.resize(0);
  return;
}

void backwardReaction::addProduct (std::string name)
{
  products.push_back(name);
  return;
}

std::string backwardReaction::getParent () const
{
  return parent;
}

std::string backwardReaction::getMonomer () const
{
  return monomer;
}

std::vector<std::string> backwardReaction::getProducts () const
{
  return products;
}
