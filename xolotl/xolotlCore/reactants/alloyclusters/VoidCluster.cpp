// Includes
#include "VoidCluster.h"
#include "AlloyClusterReactionNetwork.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

VoidCluster::VoidCluster(int n,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		AlloyCluster(registry) {
	// Set the size
	size = n;
	// Update the composition map
	compositionMap[voidType] = size;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "Vo_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = voidType;

	// Compute the reaction radius
	double FourPi = 4.0 * xolotlCore::pi;

	// DEFINE THE RADIUS HERE
	reactionRadius = 0.0;

	return;
}
