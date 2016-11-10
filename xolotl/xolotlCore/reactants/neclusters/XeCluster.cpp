// Includes
#include "XeCluster.h"
#include "NEClusterReactionNetwork.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

XeCluster::XeCluster(int nXe,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		NECluster(registry) {
	// Set the size
	size = nXe;
	// Update the composition map
	compositionMap[xeType] = size;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "Xe_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = "Xe";

	// Compute the reaction radius
	double FourPi = 4.0 * xolotlCore::pi;
	reactionRadius = 1.05 * pow((3.0 * 85.0 * (double) size) / FourPi, (1.0/3.0)) / 10.0;
	if (size == 1) reactionRadius = 0.3;

	return;
}

void XeCluster::createReactionConnectivity() {
	// Call the function from the NECluster class to take care of the single
	// species reactions
	NECluster::createReactionConnectivity();

	// This cluster is always Xe_a

	return;
}

void XeCluster::createDissociationConnectivity() {
	// Call the function from the NECluster class to take care of the single
	// species dissociation
	NECluster::createDissociationConnectivity();

	// This cluster is always Xe_a

	return;
}
