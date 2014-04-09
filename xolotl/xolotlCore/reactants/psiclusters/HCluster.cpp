// Includes
#include "HCluster.h"

using namespace xolotlCore;

HCluster::HCluster(int nH, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(nH, registry) {
	// Set the reactant name appropriately
	name = "Hydrogen";
}

HCluster::~HCluster() { }
