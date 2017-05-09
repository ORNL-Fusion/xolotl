// Includes
#include "VacCluster.h"
#include "AlloyClusterReactionNetwork.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

VacCluster::VacCluster(int n,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		AlloyCluster(registry) {
	// Set the size
	size = n;
	// Update the composition map
	compositionMap[vType] = size;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "V_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = vType;

	// Compute the reaction radius
	double FourPi = 4.0 * xolotlCore::pi;

	// DEFINE THE RADIUS HERE
	reactionRadius = 0.0;

	return;
}
