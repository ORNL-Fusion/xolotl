#include "TCluster.h"

using namespace xolotlCore;

TCluster::TCluster(int nT) :
		HCluster(nT) {
	// Set the reactant name appropriately
	name = "Tritium";
}
TCluster::~TCluster() {
}
