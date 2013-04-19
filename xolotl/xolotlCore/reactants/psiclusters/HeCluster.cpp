// Includes
#include "HeCluster.h"

using namespace xolotlCore;

HeCluster::HeCluster(int nHe) :
		PSICluster(nHe) {
	// Set the reactant name appropriately
	name = "Helium";
}
HeCluster::~HeCluster() {
}
