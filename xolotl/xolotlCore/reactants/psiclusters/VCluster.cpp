// Includes
#include "VCluster.h"

using namespace xolotlCore;

VCluster::VCluster(int nV) :
		PSICluster(nV) {
	// Set the reactant name appropriately
	name = "Vacancy";
}
VCluster::~VCluster() {
}
