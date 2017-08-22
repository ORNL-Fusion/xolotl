#include "DCluster.h"
#include <Constants.h>

using namespace xolotlCore;

DCluster::DCluster(int nD,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		HCluster(registry) {
	// Set the size
	size = nD;
	// Update the composition map
	compositionMap[dType] = size;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "D_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = dType;

	// Compute the reaction radius
	double a0 = 0.02;
	double FourPi = 4.0 * xolotlCore::pi;
	reactionRadius = a0 * pow((3.0 / FourPi) * size, (1.0 / 3.0));

	return;
}
