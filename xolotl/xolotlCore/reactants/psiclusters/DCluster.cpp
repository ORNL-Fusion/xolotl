#include "DCluster.h"

using namespace xolotlCore;

DCluster::DCluster(int nH) :
		HCluster(nH) {
	// Set the reactant name appropriately
	name = "Deuterium";
}
DCluster::~DCluster() {
}
