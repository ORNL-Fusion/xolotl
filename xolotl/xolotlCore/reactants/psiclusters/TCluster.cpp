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
	double a0 = 0.03;
	double FourPi = 4.0 * xolotlCore::pi;
	reactionRadius = a0 * pow((3.0 / FourPi) * size, (1.0 / 3.0));

	return;
}
