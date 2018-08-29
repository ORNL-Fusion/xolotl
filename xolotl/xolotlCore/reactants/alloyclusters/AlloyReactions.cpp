#include "AlloyReactions.h"

using namespace xolotlCore;

forwardReaction::forwardReaction(ReactantType name1, ReactantType name2) {
	reactant1 = name1;
	reactant2 = name2;
	products.clear();
	return;
}

void forwardReaction::addProduct(ReactantType name) {
	products.push_back(name);
	return;
}

ReactantType forwardReaction::getFirstReactant() const {
	return reactant1;
}

ReactantType forwardReaction::getSecondReactant() const {
	return reactant2;
}

std::vector<ReactantType> forwardReaction::getProducts() const {
	return products;
}

backwardReaction::backwardReaction(ReactantType parentName,
		ReactantType monomerName) {
	parent = parentName;
	monomer = monomerName;
	products.clear();
	return;
}

void backwardReaction::addProduct(ReactantType name) {
	products.push_back(name);
	return;
}

ReactantType backwardReaction::getParent() const {
	return parent;
}

ReactantType backwardReaction::getMonomer() const {
	return monomer;
}

std::vector<ReactantType> backwardReaction::getProducts() const {
	return products;
}
