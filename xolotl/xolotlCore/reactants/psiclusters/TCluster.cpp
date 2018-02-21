#include "TCluster.h"
#include <Constants.h>

using namespace xolotlCore;

TCluster::TCluster(int nT,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		HCluster(registry) {
	// Set the size
	size = nT;
	// Update the composition map
	compositionMap[tType] = size;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "T_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = tType;

	// Compute the reaction radius
	reactionRadius = 0.3;

	return;
}
