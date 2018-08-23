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

	// Define the diffusion pre-factor
	{
	  double jumpDistance = xolotlCore::alloyLatticeConstant / sqrt(2.0);
	  double phononFrequency = 1.0e13;
	  double jumpsPerPhonon = 1.0;
	  double prefactorExponent = -1.0;
	  diffusionFactor = phononFrequency * jumpsPerPhonon * jumpDistance
	      * jumpDistance * pow(double(size),prefactorExponent) / (6.0);
	}

	// Define the formation energy
	{
		double A = 1.5; // vacancy formation energy
		double B = 3.41649; // fit to Zinkle formula
		formationEnergy = A + B * ( pow(double(size),2.0/3.0) - 1.0 );
	}

	// Define the migration energy
	migrationEnergy = 1.2;

	// Define the reaction radius
	reactionRadius = 0.5 * xolotlCore::alloyLatticeConstant
	    * pow( 1.5 * double(size) / xolotlCore::pi , 1.0/3.0 );

	return;
}
