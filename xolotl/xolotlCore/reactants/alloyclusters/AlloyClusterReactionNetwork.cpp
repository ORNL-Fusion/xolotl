#include "AlloyClusterReactionNetwork.h"
#include "AlloyCluster.h"
#include "AlloySuperCluster.h"
#include <xolotlPerf.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <Constants.h>

using namespace xolotlCore;

void AlloyClusterReactionNetwork::setDefaultPropsAndNames() {
	// Shared pointers for the cluster type map
	std::shared_ptr<std::vector<std::shared_ptr<IReactant>>>vVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> iVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> voidVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> faultedVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> frankVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> perfectVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> voidSuperVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> faultedSuperVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> frankSuperVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> perfectSuperVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> alloySuperVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();

	// Initialize default properties
	dissociationsEnabled = true;
	numVClusters = 0;
	numIClusters = 0;
	numVoidClusters = 0;
	numFaultedClusters = 0;
	numFrankClusters = 0;
	numPerfectClusters = 0;
	numSuperClusters = 0;
	maxVClusterSize = 0;
	maxIClusterSize = 0;
	maxVoidClusterSize = 0;
	maxFaultedClusterSize = 0;
	maxFrankClusterSize = 0;
	maxPerfectClusterSize = 0;

	// Initialize the current and last size to 0
	networkSize = 0;
	// Set the reactant names
	names.push_back(vType);
	names.push_back(iType);
	names.push_back(voidType);
	names.push_back(faultedType);
	names.push_back(frankType);
	names.push_back(perfectType);
	names.push_back(AlloyVoidSuperType);
	names.push_back(AlloyFaultedSuperType);
	names.push_back(AlloyFrankSuperType);
	names.push_back(AlloyPerfectSuperType);
	names.push_back(AlloySuperType);

	// Setup the cluster type map
	clusterTypeMap[vType] = vVector;
	clusterTypeMap[iType] = iVector;
	clusterTypeMap[voidType] = voidVector;
	clusterTypeMap[faultedType] = faultedVector;
	clusterTypeMap[frankType] = frankVector;
	clusterTypeMap[perfectType] = perfectVector;
	clusterTypeMap[AlloyVoidSuperType] = voidSuperVector;
	clusterTypeMap[AlloyFaultedSuperType] = faultedSuperVector;
	clusterTypeMap[AlloyFrankSuperType] = frankSuperVector;
	clusterTypeMap[AlloyPerfectSuperType] = perfectSuperVector;
	clusterTypeMap[AlloySuperType] = alloySuperVector;

	return;
}

AlloyClusterReactionNetwork::AlloyClusterReactionNetwork() :
		ReactionNetwork() {
	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}

AlloyClusterReactionNetwork::AlloyClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork(registry) {
	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}

AlloyClusterReactionNetwork::AlloyClusterReactionNetwork(
		const AlloyClusterReactionNetwork &other) :
		ReactionNetwork(other) {
	// The size and ids do not need to be copied. They will be fixed when the
	// reactants are added.

	// Reset the properties table so that it can be properly updated when the
	// network is filled.
	setDefaultPropsAndNames();
	// Get all of the reactants from the other network and add them to this one
	// Load the single-species clusters. Calling getAll() will not work because
	// it is not const.
	std::vector<std::shared_ptr<IReactant> > reactants;
	for (auto it = other.singleSpeciesMap.begin();
			it != other.singleSpeciesMap.end(); ++it) {
		reactants.push_back(it->second);
	}
	// Load the mixed-species clusters
	for (auto it = other.mixedSpeciesMap.begin();
			it != other.mixedSpeciesMap.end(); ++it) {
		reactants.push_back(it->second);
	}
	// Load the super-species clusters
	for (auto it = other.superSpeciesMap.begin();
			it != other.superSpeciesMap.end(); ++it) {
		reactants.push_back(it->second);
	}
	for (unsigned int i = 0; i < reactants.size(); i++) {
		add(reactants[i]->clone());
	}

	return;
}

double AlloyClusterReactionNetwork::calculateReactionRateConstant(
		ProductionReaction * reaction) const {
	// Get the reaction radii
	double r_first = reaction->first->getReactionRadius();
	double r_second = reaction->second->getReactionRadius();

	// Get the diffusion coefficients
	double firstDiffusion = reaction->first->getDiffusionCoefficient();
	double secondDiffusion = reaction->second->getDiffusionCoefficient();

	// Calculate and return
	double k_plus = 4.0 * xolotlCore::pi
			* (r_first + r_second + xolotlCore::alloyCoreRadius)
			* (firstDiffusion + secondDiffusion);
	return k_plus;
}

double AlloyClusterReactionNetwork::calculateDissociationConstant(
		DissociationReaction * reaction) const {
	// If the dissociations are not allowed
	if (!dissociationsEnabled)
		return 0.0;

	// Compute the atomic volume (included the fact that there are 4 atoms per cell)
	double atomicVolume = 0.25 * xolotlCore::alloyLatticeConstant
			* xolotlCore::alloyLatticeConstant
			* xolotlCore::alloyLatticeConstant;

	// Get the rate constant from the reverse reaction
	double kPlus = reaction->reverseReaction->kConstant;

	// Calculate and return
	double bindingEnergy = computeBindingEnergy(reaction);
	double k_minus_exp = exp(
			-1.0 * bindingEnergy / (xolotlCore::kBoltzmann * temperature));
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;

	return k_minus;
}

void AlloyClusterReactionNetwork::createReactionConnectivity() {
	// Initial declarations
	int firstSize = 0, secondSize = 0, productSize = 0;

	// All reactions at once
	// Make a vector of types
	std::vector<string> typeVec = { vType, iType, voidType, faultedType,
			frankType, perfectType };
	// Loop on it
	for (auto tvIter1 = typeVec.begin(); tvIter1 != typeVec.end(); ++tvIter1) {
		// Get name of first reactant
		string typeName1 = *tvIter1;
		// Gett all reactants of type 1
		auto allReactants1 = getAll(typeName1);
		// Loop on it again starting at same iterator to avoid repeats
		for (auto tvIter2 = tvIter1; tvIter2 != typeVec.end(); ++tvIter2) {
			// Get name of second reactants
			string typeName2 = *tvIter2;
			// Special case with single species reaction
			if (typeName1 == typeName2) {

				// Loop on the clusters
				for (auto firstIt = allReactants1.begin();
						firstIt != allReactants1.end(); firstIt++) {
					// Get its size
					firstSize = (*firstIt)->getSize();
					// Check if type 1 is vacancy type
					if (typeName1 == vType || typeName1 == faultedType
							|| typeName1 == voidType) {
						firstSize = -firstSize;
					}
					// Loop on the clusters starting at iterator so no repeat
					for (auto secondIt = firstIt;
							secondIt != allReactants1.end(); secondIt++) {
						// Get its size
						secondSize = (*secondIt)->getSize();
						// Only proceed if at least one reactant is mobile
						if ((*firstIt)->getDiffusionFactor() > 0.0
								|| (*secondIt)->getDiffusionFactor() > 0.0) {
							// Check if type 2 is vacancy type
							if (typeName2 == vType || typeName2 == faultedType
									|| typeName2 == voidType) {
								secondSize = -secondSize;
							}
							// Find product size
							productSize = firstSize + secondSize;
							// Interstitial type product
							if (productSize > 0) {
								// Get the products
								IReactant * product = nullptr;
								auto product1 = get(iType, productSize);
								IReactant * product2 = nullptr;
								if (typeName1 == perfectType) {
									product2 = get(perfectType, productSize);
								}
								auto product3 = get(frankType, productSize);
								// Check which product if any is produced
								bool isReaction = true;
								if (product1) {
									product = get(iType, productSize);
								} else if (product2) {
									product = get(perfectType, productSize);
								} else if (product3) {
									product = get(frankType, productSize);
								} else {
									isReaction = false;
									// Throw a warning?
								}
								// Create reaction and connectivity if there is a reaction
								if (isReaction) {
									// Create a production reaction
									auto reaction = std::make_shared<
											ProductionReaction>((*firstIt),
											(*secondIt));
									// Tell the reactants that they are in this reaction
									(*firstIt)->createCombination(reaction);
									(*secondIt)->createCombination(reaction);
									product->createProduction(reaction);

									//std::cout << (*firstIt)->getName() << " + "
									//		<< (*secondIt)->getName() << " -> "
									//		<< product->getName() << std::endl;

									// Check if the reverse reaction is allowed
									checkDissociationConnectivity(product,
											reaction);
								}
							}
							// Vacancy type product
							else if (productSize < 0) {
								productSize = -productSize;
								// Get the products
								IReactant * product = nullptr;
								auto product1 = get(vType, productSize);
								IReactant * product2 = nullptr;
								if (typeName1 == faultedType) {
									product2 = get(faultedType, productSize);
								}
								auto product3 = get(voidType, productSize);
								// Check which product if any is produced
								bool isReaction = true;
								if (product1) {
									product = get(vType, productSize);
								} else if (product2) {
									product = get(faultedType, productSize);
								} else if (product3) {
									product = get(voidType, productSize);
								} else {
									isReaction = false;
									// Throw a warning?
								}
								// Create reaction and connectivity if there is a reaction
								if (isReaction) {
									// Create a production reaction
									auto reaction = std::make_shared<
											ProductionReaction>((*firstIt),
											(*secondIt));
									// Tell the reactants that they are in this reaction
									(*firstIt)->createCombination(reaction);
									(*secondIt)->createCombination(reaction);
									product->createProduction(reaction);

									//std::cout << (*firstIt)->getName() << " + "
									//		<< (*secondIt)->getName() << " -> "
									//		<< product->getName() << std::endl;

									// Check if the reverse reaction is allowed
									checkDissociationConnectivity(product,
											reaction);
								}
							}
//							// Recombination
//							else {
//								// Shouldn't happen, throw warning?
//								// Create a production reaction
//								auto reaction = std::make_shared<
//										ProductionReaction>((*firstIt),
//										(*secondIt));
//								// Tell the reactants that they are in this reaction
//								(*firstIt)->createCombination(reaction);
//								(*secondIt)->createCombination(reaction);
//
//								std::cout << (*firstIt)->getName() << " + "
//										<< (*secondIt)->getName() << std::endl;
//							}
						}
					}
				}

			}
			// All other reactions
			else {
				// Get all reactants of type 2
				auto allReactants2 = getAll(typeName2);

				// Loop on the type 1 clusters
				for (auto firstIt = allReactants1.begin();
						firstIt != allReactants1.end(); firstIt++) {
					// Get its size
					firstSize = (*firstIt)->getSize();
					// Check if type 1 is vacancy type
					if (typeName1 == vType || typeName1 == faultedType
							|| typeName1 == voidType) {
						firstSize = -firstSize;
					}
					// Loop on the type 2 clusters
					for (auto secondIt = allReactants2.begin();
							secondIt != allReactants2.end(); secondIt++) {
						// Get its size
						secondSize = (*secondIt)->getSize();
						// Only proceed if at least one reactant is mobile
						if ((*firstIt)->getDiffusionFactor() > 0.0
								|| (*secondIt)->getDiffusionFactor() > 0.0) {
							// Check if type 2 is vacancy type
							if (typeName2 == vType || typeName2 == faultedType
									|| typeName2 == voidType) {
								secondSize = -secondSize;
							}
							// Find product size
							productSize = firstSize + secondSize;
							// Interstitial type product
							if (productSize > 0) {
								// Get the products
								IReactant * product = nullptr;
								auto product1 = get(iType, productSize);
								IReactant * product2 = nullptr;
								if (typeName1 == perfectType
										|| typeName2 == perfectType) {
									product2 = get(perfectType, productSize);
								}
								auto product3 = get(frankType, productSize);
								// Check which product if any is produced
								bool isReaction = true;
								if (product1) {
									product = get(iType, productSize);
								} else if (product2) {
									product = get(perfectType, productSize);
								} else if (product3) {
									product = get(frankType, productSize);
								} else {
									isReaction = false;
									// Throw an warning?
								}
								// Create reaction and connectivity if there is a reaction
								if (isReaction) {
									// Create a production reaction
									auto reaction = std::make_shared<
											ProductionReaction>((*firstIt),
											(*secondIt));
									// Tell the reactants that they are in this reaction
									(*firstIt)->createCombination(reaction);
									(*secondIt)->createCombination(reaction);
									product->createProduction(reaction);

									//std::cout << (*firstIt)->getName() << " + "
									//		<< (*secondIt)->getName() << " -> "
									//		<< product->getName() << std::endl;

									// Check if the reverse reaction is allowed
									checkDissociationConnectivity(product,
											reaction);
								}
							}
							// Vacancy type product
							else if (productSize < 0) {
								productSize = -productSize;
								// Get the products
								IReactant * product = nullptr;
								auto product1 = get(vType, productSize);
								IReactant * product2 = nullptr;
								if ((typeName1 != voidType
										&& typeName2 != voidType)
										&& (typeName1 != vType
												|| typeName2 != vType)) {
									product2 = get(faultedType, productSize);
								}
								auto product3 = get(voidType, productSize);
								// Check which product if any is produced
								bool isReaction = true;
								if (product1) {
									product = get(vType, productSize);
								} else if (product2) {
									product = get(faultedType, productSize);
								} else if (product3) {
									product = get(voidType, productSize);
								} else {
									isReaction = false;
									// Throw a warning?
								}
								// Create reaction and connectivity if there is a reaction
								if (isReaction) {
									// Create a production reaction
									auto reaction = std::make_shared<
											ProductionReaction>((*firstIt),
											(*secondIt));
									// Tell the reactants that they are in this reaction
									(*firstIt)->createCombination(reaction);
									(*secondIt)->createCombination(reaction);
									product->createProduction(reaction);

									//std::cout << (*firstIt)->getName() << " + "
									//		<< (*secondIt)->getName() << " -> "
									//		<< product->getName() << std::endl;

									// Check if the reverse reaction is allowed
									checkDissociationConnectivity(product,
											reaction);
								}
							}
							// Recombination
							else {
								// Create a production reaction
								auto reaction = std::make_shared<
										ProductionReaction>((*firstIt),
										(*secondIt));
								// Tell the reactants that they are in this reaction
								(*firstIt)->createCombination(reaction);
								(*secondIt)->createCombination(reaction);

								//std::cout << (*firstIt)->getName() << " + "
								//		<< (*secondIt)->getName() << std::endl;
							}
						}
					}
				}

			}
		}
	}

	// One way reaction faulted loop to vacancy and vacancy cluster
	{
		auto allFaulted = getAll(faultedType);
		int minSize = 0;
		for (auto it = allFaulted.begin(); it != allFaulted.end(); ++it) {
			int faultedSize = (*it)->getSize();
			if (minSize == 0) {
				minSize = faultedSize;
			}
			else if (faultedSize < minSize) {
				minSize = faultedSize;
			}
		}
		auto product = get(faultedType, minSize);
		auto reactant1 = get(vType, 1);
		auto reactant2 = get(vType, (minSize - 1));
		if (product && reactant1 && reactant2) {
			  if (reactant1->getDiffusionFactor() > 0.0
				  	|| reactant2->getDiffusionFactor() > 0.0) {
			    auto reaction = std::make_shared<
					    ProductionReaction>((reactant1),
					    (reactant2));
			    checkDissociationConnectivity(product,
					    reaction);
			}
		}
	}


	return;
}

void AlloyClusterReactionNetwork::checkDissociationConnectivity(
		IReactant * emittingReactant,
		std::shared_ptr<ProductionReaction> reaction) {

	// Get typenames for two reactants
	string firstName = reaction->first->getType();
	string secondName = reaction->second->getType();

	// Check if at least one of the potentially emitted cluster is size one
	if (reaction->first->getSize() != 1 && reaction->second->getSize() != 1) {
		// Don't add the reverse reaction
		return;
	}
	// Only allow vacancy type reactions to be reversible
	else if (firstName == iType || firstName == frankType
			|| firstName == perfectType || secondName == iType
			|| secondName == frankType || secondName == perfectType) {
		return;
	}

	// The reaction can occur, create the dissociation
	// Create a dissociation reaction
	auto dissociationReaction = std::make_shared<DissociationReaction>(
			emittingReactant, reaction->first, reaction->second);
	// Set the reverse reaction
	dissociationReaction->reverseReaction = reaction.get();
	// Tell the reactants that their are in this reaction
	reaction->first->createDissociation(dissociationReaction);
	reaction->second->createDissociation(dissociationReaction);
	emittingReactant->createEmission(dissociationReaction);

	//std::cout << (emittingReactant)->getName() << " -> "
	//		<< reaction->first->getName() << " + "
	//		<< reaction->second->getName() << std::endl;

	return;
}

void AlloyClusterReactionNetwork::setTemperature(double temp) {
	ReactionNetwork::setTemperature(temp);

	computeRateConstants();

	return;
}

double AlloyClusterReactionNetwork::getTemperature() const {
	return temperature;
}

IReactant * AlloyClusterReactionNetwork::get(const std::string& type,
		const int size) const {
	// Local Declarations
	static std::map<std::string, int> composition = { { vType, 0 },
			{ iType, 0 }, { voidType, 0 }, { faultedType, 0 }, { frankType, 0 },
			{ perfectType, 0 } };
	std::shared_ptr<IReactant> retReactant;

	// Setup the composition map to default values because it is static
	composition[vType] = 0;
	composition[iType] = 0;
	composition[voidType] = 0;
	composition[faultedType] = 0;
	composition[frankType] = 0;
	composition[perfectType] = 0;

	// Only pull the reactant if the name and size are valid
	if ((type == vType || type == iType || type == voidType
			|| type == faultedType || type == frankType || type == perfectType)
			&& size >= 1) {
		composition[type] = size;
		//std::string encodedName = AlloyCluster::encodeCompositionAsName(composition);
		// Make sure the reactant is in the map
		std::string compStr = Reactant::toCanonicalString(type, composition);
		if (singleSpeciesMap.count(compStr)) {
			retReactant = singleSpeciesMap.at(compStr);
		}
	}

	return retReactant.get();
}

IReactant * AlloyClusterReactionNetwork::getCompound(const std::string& type,
		const std::vector<int>& sizes) const {
	// Local Declarations
	std::shared_ptr<IReactant> retReactant;

	return retReactant.get();
}

IReactant * AlloyClusterReactionNetwork::getSuper(const std::string& type,
		const int size) const {
	// Local Declarations
	static std::map<std::string, int> composition = { { vType, 0 },
			{ iType, 0 }, { voidType, 0 }, { faultedType, 0 }, { frankType, 0 },
			{ perfectType, 0 } };
	std::shared_ptr<IReactant> retReactant;

	// Setup the composition map to default values because it is static
	composition[vType] = 0;
	composition[iType] = 0;
	composition[voidType] = 0;
	composition[faultedType] = 0;
	composition[frankType] = 0;
	composition[perfectType] = 0;

	// Only pull the reactant if the name and size are valid
	if (size >= 1) {
		if (type == AlloyVoidSuperType)
			composition[voidType] = size;
		else if (type == AlloyFaultedSuperType)
			composition[faultedType] = size;
		else if (type == AlloyFrankSuperType)
			composition[frankType] = size;
		else if (type == AlloyPerfectSuperType)
			composition[perfectType] = size;
		// Make sure the reactant is in the map
		std::string compStr = Reactant::toCanonicalString(type, composition);
		if (superSpeciesMap.count(compStr)) {
			retReactant = superSpeciesMap.at(compStr);
		}
	}

	return retReactant.get();
}

const std::shared_ptr<std::vector<IReactant *>> & AlloyClusterReactionNetwork::getAll() const {
	return allReactants;
}

std::vector<IReactant *> AlloyClusterReactionNetwork::getAll(
		const std::string& name) const {
	// Local Declarations
	std::vector<IReactant *> reactants;

	// Only pull the reactants if the name is valid
	if (name == vType || name == iType || name == voidType
			|| name == faultedType || name == frankType || name == perfectType
			|| name == AlloyVoidSuperType || name == AlloyFaultedSuperType
			|| name == AlloyFrankSuperType || name == AlloyPerfectSuperType
			|| name == AlloySuperType) {
		std::shared_ptr<std::vector<std::shared_ptr<IReactant>> > storedReactants =
				clusterTypeMap.at(name);
		int vecSize = storedReactants->size();
		for (int i = 0; i < vecSize; i++) {
			reactants.push_back(storedReactants->at(i).get());
		}
	}

	return reactants;
}

void AlloyClusterReactionNetwork::add(std::shared_ptr<IReactant> reactant) {
	// Local Declarations
	int numV = 0, numI = 0, numVoid = 0, numFaulted = 0, numFrank = 0,
			numPerfect = 0;
	int* numClusters = nullptr;
	int* maxClusterSize = nullptr;

	// Only add a complete reactant
	if (reactant != NULL) {
		// Get the composition
		auto composition = reactant->getComposition();
		std::string compStr = reactant->getCompositionString();

		// Get the species sizes
		numV = composition.at(vType);
		numI = composition.at(iType);
		numVoid = composition.at(voidType);
		numFaulted = composition.at(faultedType);
		numFrank = composition.at(frankType);
		numPerfect = composition.at(perfectType);

		if (singleSpeciesMap.count(compStr) == 0) {
			/// Put the reactant in its map
			singleSpeciesMap[compStr] = reactant;

			// Figure out which type we have
			if (numV > 0) {
				numClusters = &numVClusters;
				maxClusterSize = &maxVClusterSize;
			} else if (numI > 0) {
				numClusters = &numIClusters;
				maxClusterSize = &maxIClusterSize;
			} else if (numVoid > 0) {
				numClusters = &numVoidClusters;
				maxClusterSize = &maxVoidClusterSize;
			} else if (numFaulted > 0) {
				numClusters = &numFaultedClusters;
				maxClusterSize = &maxFaultedClusterSize;
			} else if (numFrank > 0) {
				numClusters = &numFrankClusters;
				maxClusterSize = &maxFrankClusterSize;
			} else {
				numClusters = &numPerfectClusters;
				maxClusterSize = &maxPerfectClusterSize;
			}
		} else {
			std::stringstream errStream;
			errStream << "AlloyClusterReactionNetwork Message: "
					<< "Duplicate Reactant (V=" << numV << ",I=" << numI
					<< ",void=" << numVoid << ",faulted=" << numFaulted
					<< ",frank=" << numFrank << ",perfect=" << numPerfect
					<< ") not added!" << std::endl;
			throw errStream.str();
		}

		// Increment the number of total clusters of this type
		(*numClusters)++;
		// Increment the max cluster size key
		int clusterSize = numV + numI + numVoid + numFaulted + numFrank
				+ numPerfect;
		(*maxClusterSize) = std::max(clusterSize, *maxClusterSize);
		// Update the size
		++networkSize;
		// Set the id for this cluster
		reactant->setId(networkSize);
		// Get the vector for this reactant from the type map
		auto clusters = clusterTypeMap[reactant->getType()];

		clusters->push_back(reactant);
		// Add the pointer to the list of all clusters
		allReactants->push_back(reactant.get());
	}

	return;
}

void AlloyClusterReactionNetwork::addSuper(
		std::shared_ptr<IReactant> reactant) {
	// Local Declarations
	int numV = 0, numI = 0, numVoid = 0, numFaulted = 0, numFrank = 0,
			numPerfect = 0;
	int* numClusters = nullptr;

	// Only add a complete reactant
	if (reactant != NULL) {
		// Get the composition
		auto composition = reactant->getComposition();
		std::string compStr = reactant->getCompositionString();

		// Get the species sizes
		numV = composition.at(vType);
		numI = composition.at(iType);
		numVoid = composition.at(voidType);
		numFaulted = composition.at(faultedType);
		numFrank = composition.at(frankType);
		numPerfect = composition.at(perfectType);

		// Only add the element if we don't already have it
		if (superSpeciesMap.count(compStr) == 0) {
			// Put the compound in its map
			superSpeciesMap[compStr] = reactant;
			// Set the key
			numClusters = &numSuperClusters;
		} else {
			std::stringstream errStream;
			errStream << "AlloyClusterReactionNetwork Message: "
					<< "Duplicate Reactant (V=" << numV << ",I=" << numI
					<< ",void=" << numVoid << ",faulted=" << numFaulted
					<< ",frank=" << numFrank << ",perfect=" << numPerfect
					<< ") not added!" << std::endl;
			throw errStream.str();
		}

		// Increment the number of total clusters of this type
		(*numClusters)++;
		// Update the size
		++networkSize;
		// Set the id for this cluster
		reactant->setId(networkSize);
		// Get the vector for this reactant from the type map
		auto clusters = clusterTypeMap[reactant->getType()];
		clusters->push_back(reactant);
		// Get the vector for super clusters
		clusters = clusterTypeMap[AlloySuperType];
		clusters->push_back(reactant);
		// Add the pointer to the list of all clusters
		allReactants->push_back(reactant.get());
	}

	return;
}

void AlloyClusterReactionNetwork::removeReactants(
		const std::vector<IReactant*>& doomedReactants) {

	// Build a ReactantMatcher functor for the doomed reactants.
	// Doing this here allows us to construct the canonical composition
	// strings for the doomed reactants once and reuse them.
	// If we used an anonymous functor object in the std::remove_if
	// calls we would build these strings several times in this function.
	ReactionNetwork::ReactantMatcher doomedReactantMatcher(doomedReactants);

	// Remove the doomed reactants from our collection of all known reactants.
	auto ariter = std::remove_if(allReactants->begin(), allReactants->end(),
			doomedReactantMatcher);
	allReactants->erase(ariter, allReactants->end());

	// Remove the doomed reactants from the type-specific cluster vectors.
	// First, determine all cluster types used by clusters in the collection
	// of doomed reactants...
	std::set<std::string> typesUsed;
	for (auto reactant : doomedReactants) {
		typesUsed.insert(reactant->getType());
	}

	// ...Next, examine each type's collection of clusters and remove the
	// doomed reactants.
	for (auto currType : typesUsed) {
		auto clusters = clusterTypeMap[currType];
		auto citer = std::remove_if(clusters->begin(), clusters->end(),
				doomedReactantMatcher);
		clusters->erase(citer, clusters->end());
	}

	// Remove the doomed reactants from the SpeciesMap.
	// We cannot use std::remove_if and our ReactantMatcher here
	// because std::remove_if reorders the elements in the underlying
	// container to move the doomed elements to the end of the container,
	// but the std::map doesn't support reordering.
	for (auto reactant : doomedReactants) {
		if (reactant->isMixed())
			mixedSpeciesMap.erase(reactant->getCompositionString());
		else
			singleSpeciesMap.erase(reactant->getCompositionString());
	}

	return;
}

void AlloyClusterReactionNetwork::reinitializeNetwork() {
	// Reset the Ids
	int id = 0;
	for (auto it = allReactants->begin(); it != allReactants->end(); ++it) {
		id++;
		(*it)->setId(id);
		(*it)->setMomentId(id);

		(*it)->optimizeReactions();
	}

	// Reset the network size
	networkSize = id;

	// Get all the super clusters and loop on them
	for (auto it = clusterTypeMap[AlloySuperType]->begin();
			it != clusterTypeMap[AlloySuperType]->end(); ++it) {
		id++;
		(*it)->setMomentId(id);
	}

	return;
}

void AlloyClusterReactionNetwork::reinitializeConnectivities() {
	// Loop on all the reactants to reset their connectivities
	for (auto it = allReactants->begin(); it != allReactants->end(); ++it) {
		(*it)->resetConnectivities();
	}

	return;
}

void AlloyClusterReactionNetwork::updateConcentrationsFromArray(
		double * concentrations) {
	// Local Declarations
	auto reactants = getAll();
	int size = reactants->size();
	int id = 0;

	// Set the concentrations
	concUpdateCounter->increment();	// increment the update concentration counter
	for (int i = 0; i < size; i++) {
		id = reactants->at(i)->getId() - 1;
		reactants->at(i)->setConcentration(concentrations[id]);
	}

	// Set the moments
	for (int i = size - numSuperClusters; i < size; i++) {
		// Get the superCluster
		auto cluster = (AlloySuperCluster *) reactants->at(i);
		id = cluster->getId() - 1;
		cluster->setZerothMoment(concentrations[id]);
		id = cluster->getMomentId() - 1;
		cluster->setMoment(concentrations[id]);
	}

	return;
}

void AlloyClusterReactionNetwork::getDiagonalFill(int *diagFill) {
	// Get all the super clusters
	auto superClusters = getAll(AlloySuperType);

	// Degrees of freedom is the total number of clusters in the network
	const int dof = getDOF();

	// Declarations for the loop
	std::vector<int> connectivity;
	int connectivityLength, id, index;

	// Get the connectivity for each reactant
	for (int i = 0; i < networkSize; i++) {
		// Get the reactant and its connectivity
		auto reactant = allReactants->at(i);
		connectivity = reactant->getConnectivity();
		connectivityLength = connectivity.size();
		// Get the reactant id so that the connectivity can be lined up in
		// the proper column
		id = reactant->getId() - 1;
		// Create the vector that will be inserted into the dFill map
		std::vector<int> columnIds;
		// Add it to the diagonal fill block
		for (int j = 0; j < connectivityLength; j++) {
			// The id starts at j*connectivity length and is always offset
			// by the id, which denotes the exact column.
			index = id * dof + j;
			diagFill[index] = connectivity[j];
			// Add a column id if the connectivity is equal to 1.
			if (connectivity[j] == 1) {
				columnIds.push_back(j);
			}
		}
		// Update the map
		dFillMap[id] = columnIds;
	}
	// Get the connectivity for each moment
	for (int i = 0; i < superClusters.size(); i++) {
		// Get the reactant and its connectivity
		auto reactant = superClusters[i];
		connectivity = reactant->getConnectivity();
		connectivityLength = connectivity.size();
		// Get the moment id so that the connectivity can be lined up in
		// the proper column
		id = reactant->getMomentId() - 1;

		// Create the vector that will be inserted into the dFill map
		std::vector<int> columnIds;
		// Add it to the diagonal fill block
		for (int j = 0; j < connectivityLength; j++) {
			// The id starts at j*connectivity length and is always offset
			// by the id, which denotes the exact column.
			index = (id) * dof + j;
			diagFill[index] = connectivity[j];
			// Add a column id if the connectivity is equal to 1.
			if (connectivity[j] == 1) {
				columnIds.push_back(j);
			}
		}
		// Update the map
		dFillMap[id] = columnIds;
	}

	return;
}

void AlloyClusterReactionNetwork::computeRateConstants() {
	// Local declarations
	double rate = 0.0;
	// Initialize the value for the biggest production rate
	double biggestProductionRate = 0.0;

	// Loop on all the production reactions
	for (auto iter = allProductionReactions.begin();
			iter != allProductionReactions.end(); iter++) {
		// Compute the rate
		rate = calculateReactionRateConstant(iter->get());
		// Set it in the reaction
		(*iter)->kConstant = rate;

		// Check if the rate is the biggest one up to now
		if (rate > biggestProductionRate)
			biggestProductionRate = rate;
	}

	// Loop on all the dissociation reactions
	for (auto iter = allDissociationReactions.begin();
			iter != allDissociationReactions.end(); iter++) {
		// Compute the rate
		rate = calculateDissociationConstant(iter->get());
		// Set it in the reaction
		(*iter)->kConstant = rate;
	}

	// Set the biggest rate
	biggestRate = biggestProductionRate;

	return;
}

void AlloyClusterReactionNetwork::computeAllFluxes(double *updatedConcOffset) {
	// Initial declarations
	IReactant * cluster;
	AlloySuperCluster * superCluster;
	double flux = 0.0;
	int reactantIndex = 0;
	auto superClusters = getAll(AlloySuperType);

	// ----- Compute all of the new fluxes -----
	for (int i = 0; i < networkSize; i++) {
		cluster = allReactants->at(i);
		// Compute the flux
		flux = cluster->getTotalFlux();
		// Update the concentration of the cluster
		reactantIndex = cluster->getId() - 1;
		updatedConcOffset[reactantIndex] += flux;
	}

	// ---- Moments ----
	for (int i = 0; i < superClusters.size(); i++) {
		superCluster = (xolotlCore::AlloySuperCluster *) superClusters[i];

		// Compute the moment flux
		flux = superCluster->getMomentFlux();
		// Update the concentration of the cluster
		reactantIndex = superCluster->getMomentId() - 1;
		updatedConcOffset[reactantIndex] += flux;
	}

	return;
}

void AlloyClusterReactionNetwork::computeAllPartials(double *vals, int *indices,
		int *size) {
	// Initial declarations
	int reactantIndex = 0, pdColIdsVectorSize = 0;
	const int dof = getDOF();
	std::vector<double> clusterPartials;
	clusterPartials.resize(dof, 0.0);
	// Get the super clusters
	auto superClusters = getAll(AlloySuperType);

	// Update the column in the Jacobian that represents each normal reactant
//	for (int i = 0; i < networkSize - superClusters.size(); i++) {
	for (int i = 0; i < networkSize; i++) {
		auto reactant = allReactants->at(i);
		// Get the reactant index
		reactantIndex = reactant->getId() - 1;

		// Get the partial derivatives
		reactant->getPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		auto pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];

			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}
	}

	// Update the column in the Jacobian that represents the moment for the super clusters
	for (int i = 0; i < superClusters.size(); i++) {
		auto reactant = (AlloySuperCluster *) superClusters[i];

		// Get the super cluster index
		reactantIndex = reactant->getId() - 1;

		// Get the partial derivatives
		reactant->getPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		auto pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];
			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}

		// Get the moment index
		reactantIndex = reactant->getMomentId() - 1;

		// Get the partial derivatives
		reactant->getMomentPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];
			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}
	}

	return;
}
