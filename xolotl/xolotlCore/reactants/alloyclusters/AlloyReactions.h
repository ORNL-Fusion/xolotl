#ifndef ALLOYREACTIONS_H
#define ALLOYREACTIONS_H

// Includes
#include <vector>
#include <string>

namespace xolotlCore {

class forwardReaction {
private:
  std::string
    reactant1,
    reactant2;
  std::vector<std::string>
    products;
public:
  forwardReaction (std::string, std::string);
  void addProduct (std::string);
  std::string getFirstReactant () const;
  std::string getSecondReactant () const;
  std::vector<std::string> getProducts () const;
};

class backwardReaction {
private:
  std::string
    parent,
    monomer;
  std::vector<std::string>
    products;
public:
  backwardReaction (std::string, std::string);
  void addProduct (std::string);
  std::string getParent () const;
  std::string getMonomer () const;
  std::vector<std::string> getProducts () const;
};

} /* end namespace xolotlCore */
#endif
