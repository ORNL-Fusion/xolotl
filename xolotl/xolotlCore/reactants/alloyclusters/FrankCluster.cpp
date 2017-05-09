// Includes
#include "FrankCluster.h"
#include "AlloyClusterReactionNetwork.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

FrankCluster::FrankCluster(int n,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		AlloyCluster(registry) {
	// Set the size
	size = n;
	// Update the composition map
	compositionMap[frankType] = size;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "Fr_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = frankType;

	// Compute the reaction radius
	double FourPi = 4.0 * xolotlCore::pi;

	// DEFINE THE RADIUS HERE
	reactionRadius = 0.0;

	return;
}
