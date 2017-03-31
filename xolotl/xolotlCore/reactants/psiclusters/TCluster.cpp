#include "TCluster.h"

using namespace xolotlCore;

TCluster::TCluster(int nT, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		HCluster(nT, registry) {
	// Set the reactant name appropriately
	name = "Tritium";
}
TCluster::~TCluster() {
}
