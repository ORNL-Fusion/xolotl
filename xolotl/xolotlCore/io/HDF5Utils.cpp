#include "HDF5Utils.h"
#include <iostream>
#include <sstream>
#include <hdf5.h>

using namespace xolotlCore;

hid_t plistId, fileId, concGroupId, subConcGroupId, concSId, networkGroupId,
		networkSId, headerGroupId;
herr_t status;

void HDF5Utils::initializeFile(std::string fileName, int networkSize,
		int gridSize) {
	// Set up file access property list with parallel I/O access
	plistId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Create the file
	fileId = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistId);

	// Close the property list
	status = H5Pclose(plistId);

	// Create the group where the header will be stored
	headerGroupId = H5Gcreate2(fileId, "headerGroup", H5P_DEFAULT, H5P_DEFAULT,
	H5P_DEFAULT);

	// Create the group where the network will be stored
	networkGroupId = H5Gcreate2(fileId, "networkGroup", H5P_DEFAULT,
	H5P_DEFAULT, H5P_DEFAULT);

	// Create the dataspace for the network with dimension dims
	hsize_t dims[2];
	dims[0] = networkSize;
	dims[1] = 8;
	networkSId = H5Screate_simple(2, dims, NULL);

	// Create the group where the concentrations will be stored
	concGroupId = H5Gcreate2(fileId, "concentrationsGroup", H5P_DEFAULT,
	H5P_DEFAULT, H5P_DEFAULT);

	// Create, write, and close the last written time step attribute
	int lastTimeStep = -1;
	hid_t lastSId = H5Screate(H5S_SCALAR);
	hid_t lastAId = H5Acreate2(concGroupId, "lastTimeStep", H5T_STD_I32LE,
			lastSId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(lastAId, H5T_STD_I32LE, &lastTimeStep);
	status = H5Aclose(lastAId);

	return;
}

void HDF5Utils::openFile(std::string fileName) {
	// Set up file access property list with parallel I/O access
	plistId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, plistId);

	// Close the property list
	status = H5Pclose(plistId);

	// Open the concentration group
	concGroupId = H5Gopen(fileId, "concentrationsGroup", H5P_DEFAULT);
	return;
}

void HDF5Utils::fillHeader(int physicalDim, int refinement) {
	// Create, write, and close the physicalDim attribute
	hid_t dimSId = H5Screate(H5S_SCALAR);
	hid_t dimAId = H5Acreate2(headerGroupId, "physicalDim", H5T_STD_I32LE,
			dimSId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(dimAId, H5T_STD_I32LE, &physicalDim);
	status = H5Aclose(dimAId);

	// Create, write, and close the refinement attribute
	hid_t refineSId = H5Screate(H5S_SCALAR);
	hid_t refineAId = H5Acreate2(headerGroupId, "refinement", H5T_STD_I32LE,
			refineSId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(refineAId, H5T_STD_I32LE, &refinement);
	status = H5Aclose(refineAId);

	return;
}

void HDF5Utils::fillNetwork(
		std::shared_ptr<PSIClusterReactionNetwork> network) {
	// Create the array that will store the network
	int networkSize = network->size();
	double networkArray[networkSize][8];

	// Get all the reactants
	auto reactants = network->getAll();

	// Loop on them
	for (int i = 0; i < networkSize; i++) {

		// Get the i-th reactant
		auto reactant = (PSICluster *) reactants->at(i);

		// Get its composition to store it
		auto composition = reactant->getComposition();
		networkArray[i][0] = composition["He"];
		networkArray[i][1] = composition["V"];
		networkArray[i][2] = composition["I"];

		// Get its binding energies to store them
		auto bindingEnergies = reactant->getBindingEnergies();
		networkArray[i][3] = bindingEnergies.at(0); // Helium binding energy
		networkArray[i][4] = bindingEnergies.at(1); // Vacancy binding energy
		networkArray[i][5] = bindingEnergies.at(2); // Interstitial binding energy

		// Get its migration energy to store it
		double migrationEnergy = reactant->getMigrationEnergy();
		networkArray[i][6] = migrationEnergy;

		// Get its diffusion factor to store it
		double diffusionFactor = reactant->getDiffusionFactor();
		networkArray[i][7] = diffusionFactor;
	}

	// Create the dataset for the network
	hid_t datasetId = H5Dcreate2(networkGroupId, "network", H5T_IEEE_F64LE,
			networkSId,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write networkArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &networkArray);

	// Create the attribute for the network size
	hid_t networkSizeSId = H5Screate(H5S_SCALAR);
	hid_t networkSizeAId = H5Acreate2(datasetId, "networkSize", H5T_STD_I32LE,
			networkSizeSId,
			H5P_DEFAULT, H5P_DEFAULT);

	// Write it
	status = H5Awrite(networkSizeAId, H5T_STD_I32LE, &networkSize);

	// Close everything
	status = H5Sclose(networkSId);
	status = H5Sclose(networkSizeSId);
	status = H5Aclose(networkSizeAId);
	status = H5Dclose(datasetId);

	return;
}

void HDF5Utils::addConcentrationSubGroup(int timeStep, int networkSize,
		double time, double deltaTime) {
	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentration_" << timeStep;

	// Create the subgroup where the concentrations at this time step will be stored
	subConcGroupId = H5Gcreate2(concGroupId, subGroupName.str().c_str(),
			H5P_DEFAULT,
			H5P_DEFAULT, H5P_DEFAULT);

	// Create, write, and close the absolute time attribute
	hid_t timeSId = H5Screate(H5S_SCALAR);
	hid_t timeAId = H5Acreate2(subConcGroupId, "absoluteTime", H5T_IEEE_F64LE,
			timeSId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(timeAId, H5T_IEEE_F64LE, &time);
	status = H5Sclose(timeSId);
	status = H5Aclose(timeAId);

	// Create, write, and close the timestep time attribute
	hid_t deltaSId = H5Screate(H5S_SCALAR);
	hid_t deltaAId = H5Acreate2(subConcGroupId, "deltaTime", H5T_IEEE_F64LE,
			deltaSId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(deltaAId, H5T_IEEE_F64LE, &deltaTime);
	status = H5Sclose(deltaSId);
	status = H5Aclose(deltaAId);

	// Overwrite the last time step attribute of the concentration group
	hid_t lastSId = H5Screate(H5S_SCALAR);
	hid_t lastAId = H5Aopen(concGroupId, "lastTimeStep", H5P_DEFAULT);
	status = H5Awrite(lastAId, H5T_STD_I32LE, &timeStep);
	status = H5Sclose(lastSId);
	status = H5Aclose(lastAId);

	// Create property list for independent dataset write
	// (needed to be able to write the datasets without having
	// HDF5 screaming).
	plistId = H5Pcreate(H5P_DATASET_XFER);
	status = H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

	return;
}

void HDF5Utils::addConcentrationDataset(int index, int size) {
	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "position_" << index;

	// Create the dataspace for the dataset with dimension dims
	hsize_t dims[2];
	dims[0] = size;
	dims[1] = 2;
	concSId = H5Screate_simple(2, dims, NULL);

	// Create the dataset of concentrations for this position
	hid_t datasetId = H5Dcreate2(subConcGroupId, datasetName.str().c_str(),
	H5T_IEEE_F64LE, concSId,
	H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Create and write the attribute for length of the dataset
	hid_t lenSId = H5Screate(H5S_SCALAR);
	hid_t lenAId = H5Acreate2(datasetId, "datasetLength", H5T_STD_I32LE, lenSId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(lenAId, H5T_STD_I32LE, &size);

	// Create property list for independent dataset write.
	plistId = H5Pcreate(H5P_DATASET_XFER);
	status = H5Pset_dxpl_mpio(plistId, H5FD_MPIO_INDEPENDENT);

	// Close everything
	status = H5Sclose(lenSId);
	status = H5Aclose(lenAId);
	status = H5Sclose(concSId);
	status = H5Dclose(datasetId);

	return;
}

void HDF5Utils::fillConcentrations(std::vector<std::vector<double> > concVector,
		int index) {
	// Create the concentration array
	double concArray[concVector.size()][2];

	// Fill it with the concentration vector
	for (int i = 0; i < concVector.size(); i++) {
		concArray[i][0] = concVector.at(i).at(0);
		concArray[i][1] = concVector.at(i).at(1);
	}

	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "position_" << index;

	// Open the already created dataset of concentrations for this position
	hid_t datasetId = H5Dopen(subConcGroupId, datasetName.str().c_str(),
			H5P_DEFAULT);

	// Write concArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, plistId,
			&concArray);

	// Close everything
	status = H5Dclose(datasetId);

	return;
}

void HDF5Utils::finalizeFile() {
	// Close everything
	status = H5Gclose(headerGroupId);
	status = H5Gclose(networkGroupId);
	status = H5Gclose(concGroupId);
	status = H5Fclose(fileId);

	return;
}

void HDF5Utils::closeFile() {
	// Close everything
	status = H5Pclose(plistId);
	status = H5Gclose(subConcGroupId);
	status = H5Gclose(concGroupId);
	status = H5Fclose(fileId);

	return;
}

void HDF5Utils::readHeader(std::string fileName, int & physicalDim) {
	// Set up file access property list with parallel I/O access
	plistId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, plistId);

	// Close the property list
	status = H5Pclose(plistId);

	// Open the header group
	hid_t groupId = H5Gopen(fileId, "/headerGroup", H5P_DEFAULT);

	// Open and read the physicalDim attribute
	hid_t physicalDimAId = H5Aopen(groupId, "physicalDim", H5P_DEFAULT);
	status = H5Aread(physicalDimAId, H5T_STD_I32LE, &physicalDim);
	status = H5Aclose(physicalDimAId);

	// Close everything
	status = H5Gclose(groupId);
	status = H5Fclose(fileId);

	return;
}

bool HDF5Utils::hasConcentrationGroup(std::string fileName,
		int & lastTimeStep) {
	// Initialize the boolean to return
	bool hasGroup = true;

	// Set up file access property list with parallel I/O access
	plistId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, plistId);

	// Close the property list
	status = H5Pclose(plistId);

	// Check the group
	bool groupExist = H5Lexists(fileId, "/concentrationsGroup", H5P_DEFAULT);
	// If the group exist
	if (groupExist) {
		// Open the concentration group
		concGroupId = H5Gopen(fileId, "/concentrationsGroup", H5P_DEFAULT);

		// Open and read the lastTimeStep attribute
		hid_t lastAId = H5Aopen(concGroupId, "lastTimeStep", H5P_DEFAULT);
		status = H5Aread(lastAId, H5T_STD_I32LE, &lastTimeStep);
		status = H5Aclose(lastAId);

		status = H5Gclose(concGroupId);

		// if lastTimeStep is still negative the group is not valid
		if (lastTimeStep < 0)
			hasGroup = false;
	}
	// if not
	else
		hasGroup = false;

	// Close everything
	status = H5Fclose(fileId);

	return hasGroup;
}

void HDF5Utils::readTimes(std::string fileName, int lastTimeStep, double & time,
		double & deltaTime) {
	// Set up file access property list with parallel I/O access
	plistId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, plistId);

	// Close the property list
	status = H5Pclose(plistId);

	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentrationsGroup/concentration_" << lastTimeStep;

	// Open this specific concentration sub group
	subConcGroupId = H5Gopen(fileId, subGroupName.str().c_str(), H5P_DEFAULT);

	// Open and read the absoluteTime attribute
	hid_t timeAId = H5Aopen(subConcGroupId, "absoluteTime", H5P_DEFAULT);
	status = H5Aread(timeAId, H5T_IEEE_F64LE, &time);
	status = H5Aclose(timeAId);

	// Open and read the deltaTime attribute
	hid_t deltaTimeAId = H5Aopen(subConcGroupId, "deltaTime", H5P_DEFAULT);
	status = H5Aread(deltaTimeAId, H5T_IEEE_F64LE, &deltaTime);
	status = H5Aclose(deltaTimeAId);

	// Close everything
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return;
}

std::vector<std::vector<double> > HDF5Utils::readNetwork(std::string fileName) {
	// Set up file access property list with parallel I/O access
	plistId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, plistId);

	// Close the property list
	status = H5Pclose(plistId);

	// Open the dataset
	hid_t datasetId = H5Dopen(fileId, "/networkGroup/network", H5P_DEFAULT);

	// Open and read the networkSize attribute
	hid_t networkSizeAId = H5Aopen(datasetId, "networkSize", H5P_DEFAULT);
	int networkSize = 0;
	status = H5Aread(networkSizeAId, H5T_STD_I32LE, &networkSize);
	status = H5Aclose(networkSizeAId);

	// Create the array that will receive the network
	double networkArray[networkSize][8];

	// Read the data set
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			&networkArray);

	// Fill the vector to return with the dataset
	std::vector<std::vector<double> > networkVector;
	// Loop on the size of the network
	for (int i = 0; i < networkSize; i++) {
		// Create the line to give to the vector
		std::vector<double> line;
		for (int j = 0; j < 8; j++) {
			line.push_back(networkArray[i][j]);
		}
		networkVector.push_back(line);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Fclose(fileId);

	return networkVector;
}

std::vector< std::vector<double> > HDF5Utils::readGridPoint(std::string fileName,
		int lastTimeStep, int i) {
	// Create te vector to return
	std::vector< std::vector<double> > toReturn;

	// Set up file access property list with parallel I/O access
	plistId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plistId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, plistId);

	// Close the property list
	status = H5Pclose(plistId);

	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "concentrationsGroup/concentration_" << lastTimeStep
			<< "/position_" << i;

	// Check the dataset
	bool datasetExist = H5Lexists(fileId, datasetName.str().c_str(),
			H5P_DEFAULT);
	// If the dataset exists
	if (datasetExist) {
		// Open the dataset
		hid_t datasetId = H5Dopen(fileId, datasetName.str().c_str(),
				H5P_DEFAULT);

		// Read the dataset length attribute
		int length = -1;
		hid_t lenAId = H5Aopen(datasetId, "datasetLength", H5P_DEFAULT);
		status = H5Aread(lenAId, H5T_STD_I32LE, &length);
		status = H5Aclose(lenAId);

		// Create the array that will receive the concentrations
		double conc[length][2];

		// Read the data set
		status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
				H5P_DEFAULT, &conc);

		// Loop on the length
		for (int i = 0; i < length; i++) {
			// Create the concentration vector for this cluster
			std::vector<double> tmp;
			tmp.push_back(conc[i][0]);
			tmp.push_back(conc[i][1]);

			// Add it to the main vector
			toReturn.push_back(tmp);
		}

		// Close everything
		status = H5Dclose(datasetId);
	}

	status = H5Fclose(fileId);

	return toReturn;
}
