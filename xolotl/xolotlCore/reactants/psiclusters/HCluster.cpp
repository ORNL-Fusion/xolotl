// Includes
#include "HCluster.h"

using namespace xolotlCore;

HCluster::HCluster(int nH) : PSICluster(nH) {
	// Set the reactant name appropriately
	name = "Hydrogen";
}
HCluster::~HCluster() {}
