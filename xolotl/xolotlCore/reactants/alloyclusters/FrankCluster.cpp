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

	// Define the diffusion pre-factor
	{
	  //double jumpDistance = xolotlCore::alloyLatticeConstant / sqrt(2.0);
	  //double phononFrequency = 1.0e13;
	  //double jumpsPerPhonon = 1.0;
	  //double prefactorExponent = -1.0;
	  //diffusionFactor = phononFrequency * jumpsPerPhonon * jumpDistance
	  //    * jumpDistance * pow(double(size),prefactorExponent) / (6.0);
		diffusionFactor = 0.0;
	}

	// Define the formation energy
	{
		// USING PLACEHOLDER VALUES => UPDATE WITH CORRECT VALUES
		double A = 4.0;
		double B = 2.0;
		formationEnergy = A + B * ( pow(double(size),2.0/3.0) - 1.0 );
	}

	// Define the migration energy
	migrationEnergy = 0.5;

	// Define the reaction radius
	reactionRadius = 0.5 * xolotlCore::alloyLatticeConstant
	    * sqrt( double(size) * sqrt(3.0) / xolotlCore::pi );

	return;
}
