// Includes
#include "HeVCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>

using namespace xolotlCore;

HeVCluster::HeVCluster(int numT, int numD, int numHe, int numV,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(registry), numT(numT), numD(numD), numHe(numHe), numV(numV) {
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = numT + numD + numHe + numV;

	// Update the composition map
	compositionMap[tType] = numT;
	compositionMap[dType] = numD;
	compositionMap[heType] = numHe;
	compositionMap[vType] = numV;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	if (numT > 0) nameStream << "T_" << numT;
	if (numD > 0) nameStream << "D_" << numD;
	if (numHe > 0) nameStream << "He_" << numHe;
	nameStream << "V_" << numV;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = heVType;

	// Compute the reaction radius
	reactionRadius = (sqrt(3.0) / 4.0) * xolotlCore::tungstenLatticeConstant
			+ pow(
					(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0) * numV)
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0))
			- pow(
					(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0))
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0));

	return;
}

HeVCluster::HeVCluster(HeVCluster &other) :
		PSICluster(other) {
	numT = other.numT;
	numD = other.numD;
	numHe = other.numHe;
	numV = other.numV;

	return;
}
