#include "DCluster.h"

using namespace xolotlCore;

DCluster::DCluster(int nH, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		HCluster(nH, registry) {
	// Set the reactant name appropriately
	name = "Deuterium";
}
DCluster::~DCluster() {
}
