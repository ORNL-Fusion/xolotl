#ifndef ALLOYREACTIONS_H
#define ALLOYREACTIONS_H

// Includes
#include <vector>
#include <string>
#include <ReactantType.h>

namespace xolotlCore {

class forwardReaction {
private:
	ReactantType reactant1, reactant2;
	std::vector<ReactantType> products;
public:
	forwardReaction(ReactantType, ReactantType);
	void addProduct(ReactantType);
	ReactantType getFirstReactant() const;
	ReactantType getSecondReactant() const;
	std::vector<ReactantType> getProducts() const;
};

class backwardReaction {
private:
	ReactantType parent, monomer;
	std::vector<ReactantType> products;
public:
	backwardReaction(ReactantType, ReactantType);
	void addProduct(ReactantType);
	ReactantType getParent() const;
	ReactantType getMonomer() const;
	std::vector<ReactantType> getProducts() const;
};

} /* end namespace xolotlCore */
#endif
