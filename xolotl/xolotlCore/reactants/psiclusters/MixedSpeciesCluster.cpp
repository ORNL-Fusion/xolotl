// Includes
#include "MixedSpeciesCluster.h"

using namespace xolotlCore;

MixedSpeciesCluster::MixedSpeciesCluster(std::map<std::string,int> speciesMap) : PSICluster(1)
{
	numD = 0;
	numH = 0;
	numHe = 0;
	numInt = 0;
	numT = 0;
	numVac = 0;
    size = 5;
	// Set the reactant name appropriately
	name = "Mixed Species";
}
MixedSpeciesCluster::~MixedSpeciesCluster()
{
    //TODO Auto-generated method stub
}
double MixedSpeciesCluster::getGenByEm() 
{
    //TODO Auto-generated method stub
    return 0;
}
double MixedSpeciesCluster::getAnnByEm() 
{
    //TODO Auto-generated method stub
    return 0;
}
int MixedSpeciesCluster::getSpeciesSize(const std::string speciesName)
{
    //TODO Auto-generated method stub
    return 0;
}
