// Includes
#include "HeCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

HeCluster::HeCluster(int nHe,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(nHe, registry) {
	// Update the composition map
	compositionMap["He"] = size;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = "He";

	// Compute the reaction radius
	double FourPi = 4.0 * xolotlCore::pi;
	double aCubed = pow(xolotlCore::latticeConstant, 3);
	double termOne = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * size,
			(1.0 / 3.0));
	double termTwo = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed, (1.0 / 3.0));
	reactionRadius = 0.3 + termOne - termTwo;

	return;
}

HeCluster::~HeCluster() {
}

std::shared_ptr<Reactant> HeCluster::clone() {
	std::shared_ptr<Reactant> reactant(new HeCluster(*this));

	return reactant;
}

void HeCluster::combineClusters(std::vector<Reactant *> & clusters,
		std::string productName) {
	// Initial declarations
	std::map<std::string, int> secondComposition;

	// Loop on the potential combining reactants
	for (int i = 0; i < clusters.size(); i++) {
		// Get the second reactant, its composition and its index
		auto secondCluster = (PSICluster *) clusters[i];
		secondComposition = secondCluster->getComposition();
		// Check that the simple product [He_(a+c)](V_b) doesn't exist
		// b can be 0 so the simple product would be a helium cluster
		PSICluster * simpleProduct;
		if (secondComposition[vType] == 0) {
			simpleProduct = (PSICluster *) network->get(heType, size + secondComposition[heType]);
		}
		else {
			std::vector<int> comp = {size + secondComposition[heType],
				secondComposition[vType],
				secondComposition[iType]};
			simpleProduct = (PSICluster *) network->getCompound(productName, comp);
		}
		if (simpleProduct) continue;
		// The simple product doesn't exist so it will go though trap-mutation
		// The reaction is
		// (He_c)(V_b) + He_a --> [He_(a+c)][V_(b+1)] + I
		std::vector<int> comp = {size + secondComposition[heType],
				secondComposition[vType] + 1,
				secondComposition[iType]};
		auto firstProduct = (PSICluster *) network->getCompound(productName, comp);
		auto secondProduct = (PSICluster *) network->get(iType, 1);
		// If both products exist
		if (firstProduct && secondProduct) {
			// This cluster combines with the second reactant
			setReactionConnectivity(secondCluster->getId());
			// Creates the combining cluster
			// The reaction constant will be computed later and is set to 0.0 for now
			CombiningCluster combCluster(secondCluster, 0.0);
			// Push the product onto the list of clusters that combine with this one
			combiningReactants.push_back(combCluster);
		}

		// Case with I_2
		// (He_c)(V_b) + He_a --> [He_(a+c)][V_(b+2)] + I_2
		// If [He_(a+c)][V_(b+1)] does not exist
		if (firstProduct) continue;

		// Get the new products [He_(a+c)][V_(b+2)] and I_2
		comp = {size + secondComposition[heType],
				secondComposition[vType] + 2,
				secondComposition[iType]};
		firstProduct = (PSICluster *) network->getCompound(productName, comp);
		secondProduct = (PSICluster *) network->get(iType, 2);
		// If both products exist
		if (firstProduct && secondProduct) {
			// This cluster combines with the second reactant
			setReactionConnectivity(secondCluster->getId());
			// Creates the combining cluster
			// The reaction constant will be computed later and is set to 0.0 for now
			CombiningCluster combCluster(secondCluster, 0.0);
			// Push the product onto the list of clusters that combine with this one
			combiningReactants.push_back(combCluster);
		}
	}

	return;
}

void HeCluster::createReactionConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species reactions
	PSICluster::createReactionConnectivity();

	// This cluster is always He_a

	// Helium-Vacancy clustering
	// He_a + V_b --> (He_a)(V_b)
	// Get all the V clusters from the network
	auto reactants = network->getAll(vType);
	// combineClusters handles V combining with He to form HeV
	PSICluster::combineClusters(reactants, heVType);

	// Helium-Interstitial clustering
	// He_a + I_b --> (He_a)(I_b)
	// Get all the I clusters from the network
	reactants = network->getAll(iType);
	// combineClusters handles I combining with He to form HeI
	PSICluster::combineClusters(reactants, heIType);

	// Helium absorption by HeV clusters
	// He_a + (He_b)(V_c) --> [He_(a+b)](V_c)
	// Get all the HeV clusters from the network
	reactants = network->getAll(heVType);
	// combineClusters handles HeV combining with He to form HeV
	PSICluster::combineClusters(reactants, heVType);

	// Helium absorption by HeI clusters
	// He_a + (He_b)(I_c) --> [He_(a+b)](I_c)
	// Get all the HeI clusters from the network
	reactants = network->getAll(heIType);
	// combineClusters handles HeI combining with He to form HeI
	PSICluster::combineClusters(reactants, heIType);

	// Helium absorption leading to trap mutation
	// (He_c)(V_b) + He_a --> [He_(a+c)][V_(b+1)] + I
	// or
	// (He_c)(V_b) + He_a --> [He_(a+c)][V_(b+2)] + I_2
	// Get all the HeV clusters from the network
	reactants = network->getAll(heVType);
	// HeCluster::combineClusters handles He combining with HeV to go through trap-mutation
	combineClusters(reactants, heVType);
	// b can be 0 so He clusters can combine with He clusters leading to trap-mutation
	// Get all the He clusters from the network
	reactants = network->getAll(heType);
	combineClusters(reactants, heVType);

	return;
}

void HeCluster::createDissociationConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species dissociation
	PSICluster::createDissociationConnectivity();

	// This cluster is always He_a

	// Specific case for the single species cluster
	if (size == 1) {
		// He dissociation of HeV cluster is handled here
		// (He_b)(V_c) --> [He_(b-a)](V_c) + He_a
		// for a = 1
		// Get all the HeV clusters of the network
		auto allHeVReactants = network->getAll(heVType);
		for (int i = 0; i < allHeVReactants.size(); i++) {
			auto cluster = (PSICluster *) allHeVReactants[i];

			// (He_b)(V_c) is the dissociating one, [He_(b-a)](V_c) is the one
			// that is also emitted during the dissociation
			auto comp = cluster->getComposition();
			std::vector<int> compositionVec = { comp[heType] - 1, comp[vType],
					0 };
			auto smallerReactant = (PSICluster *) network->getCompound(heVType, compositionVec);
			dissociateCluster(cluster, smallerReactant);
		}

		// He dissociation of HeI cluster is handled here
		// (He_b)(I_c) --> [He_(b-a)](I_c) + He_a
		// for a = 1
		// Get all the HeI clusters of the network
		auto allHeIReactants = network->getAll(heIType);
		for (int i = 0; i < allHeIReactants.size(); i++) {
			auto cluster = (PSICluster *) allHeIReactants[i];

			// (He_b)(I_c) is the dissociating one, [He_(b-a)](I_c) is the one
			// that is also emitted during the dissociation
			auto comp = cluster->getComposition();
			std::vector<int> compositionVec = { comp[heType] - 1, 0,
					comp[iType] };
			auto smallerReactant = (PSICluster *) network->getCompound(heIType, compositionVec);
			dissociateCluster(cluster, smallerReactant);
		}
	}

	return;
}
