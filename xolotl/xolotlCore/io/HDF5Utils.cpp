#include "HDF5Utils.h"
#include <iostream>
#include <sstream>
#include <hdf5.h>
#include <mpi.h>

using namespace xolotlCore;

hid_t propertyListId, fileId, concentrationGroupId, subConcGroupId, concDataspaceId, networkGroupId,
		networkDataspaceId, headerGroupId;
herr_t status;

void HDF5Utils::initializeFile(std::string fileName, int networkSize) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Create the file
	fileId = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Create the group where the header will be stored
	headerGroupId = H5Gcreate2(fileId, "headerGroup", H5P_DEFAULT, H5P_DEFAULT,
	H5P_DEFAULT);

	// Create the group where the network will be stored
	networkGroupId = H5Gcreate2(fileId, "networkGroup", H5P_DEFAULT,
	H5P_DEFAULT, H5P_DEFAULT);

	// Create the dataspace for the network with dimension dims
	hsize_t dims[2];
	dims[0] = networkSize;
	dims[1] = 6;
	networkDataspaceId = H5Screate_simple(2, dims, NULL);

	// Create the group where the concentrations will be stored
	concentrationGroupId = H5Gcreate2(fileId, "concentrationsGroup", H5P_DEFAULT,
	H5P_DEFAULT, H5P_DEFAULT);

	// Create, write, and close the last written time step attribute
	int lastTimeStep = -1;
	hid_t lastDataspaceId = H5Screate(H5S_SCALAR);
	hid_t lastAttributeId = H5Acreate2(concentrationGroupId, "lastTimeStep", H5T_STD_I32LE,
			lastDataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(lastAttributeId, H5T_STD_I32LE, &lastTimeStep);
	status = H5Aclose(lastAttributeId);

	return;
}

void HDF5Utils::openFile(std::string fileName) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Open the concentration group
	concentrationGroupId = H5Gopen(fileId, "concentrationsGroup", H5P_DEFAULT);
	return;
}

void HDF5Utils::fillHeader(int nx, double hx, int ny,
		double hy, int nz, double hz) {
	// Create, write, and close the nx attribute
	hid_t dataspaceId = H5Screate(H5S_SCALAR);
	hid_t attributeId = H5Acreate2(headerGroupId, "nx", H5T_STD_I32LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &nx);
	status = H5Aclose(attributeId);
	// Create, write, and close the hx attribute
	attributeId = H5Acreate2(headerGroupId, "hx", H5T_IEEE_F64LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &hx);
	status = H5Aclose(attributeId);

	// Create, write, and close the ny attribute
	attributeId = H5Acreate2(headerGroupId, "ny", H5T_STD_I32LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &ny);
	status = H5Aclose(attributeId);
	// Create, write, and close the hy attribute
	attributeId = H5Acreate2(headerGroupId, "hy", H5T_IEEE_F64LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &hy);
	status = H5Aclose(attributeId);

	// Create, write, and close the nz attribute
	attributeId = H5Acreate2(headerGroupId, "nz", H5T_STD_I32LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &nz);
	status = H5Aclose(attributeId);
	// Create, write, and close the hz attribute
	attributeId = H5Acreate2(headerGroupId, "hz", H5T_IEEE_F64LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &hz);
	status = H5Aclose(attributeId);

	return;
}

void HDF5Utils::fillNetwork(PSIClusterReactionNetwork *network) {
	// Create the array that will store the network
	int networkSize = network->size();
	double networkArray[networkSize][6];

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

		// Get its formation energy to store it
		auto formationEnergy = reactant->getFormationEnergy();
		networkArray[i][3] = formationEnergy;

		// Get its migration energy to store it
		double migrationEnergy = reactant->getMigrationEnergy();
		networkArray[i][4] = migrationEnergy;

		// Get its diffusion factor to store it
		double diffusionFactor = reactant->getDiffusionFactor();
		networkArray[i][5] = diffusionFactor;
	}

	// Create the dataset for the network
	hid_t datasetId = H5Dcreate2(networkGroupId, "network", H5T_IEEE_F64LE,
			networkDataspaceId,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write networkArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &networkArray);

	// Create the attribute for the network size
	hid_t networkSizeDataspaceId = H5Screate(H5S_SCALAR);
	hid_t networkSizeAttributeId = H5Acreate2(datasetId, "networkSize", H5T_STD_I32LE,
			networkSizeDataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);

	// Write it
	status = H5Awrite(networkSizeAttributeId, H5T_STD_I32LE, &networkSize);

	// Close everything
	status = H5Sclose(networkDataspaceId);
	status = H5Sclose(networkSizeDataspaceId);
	status = H5Aclose(networkSizeAttributeId);
	status = H5Dclose(datasetId);

	return;
}

void HDF5Utils::addConcentrationSubGroup(int timeStep, int networkSize,
		double time, double deltaTime) {
	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentration_" << timeStep;

	// Create the subgroup where the concentrations at this time step will be stored
	subConcGroupId = H5Gcreate2(concentrationGroupId, subGroupName.str().c_str(),
			H5P_DEFAULT,
			H5P_DEFAULT, H5P_DEFAULT);

	// Create, write, and close the absolute time attribute
	hid_t timeDataspaceId = H5Screate(H5S_SCALAR);
	hid_t timeAttributeId = H5Acreate2(subConcGroupId, "absoluteTime", H5T_IEEE_F64LE,
			timeDataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(timeAttributeId, H5T_IEEE_F64LE, &time);
	status = H5Sclose(timeDataspaceId);
	status = H5Aclose(timeAttributeId);

	// Create, write, and close the timestep time attribute
	hid_t deltaDataspaceId = H5Screate(H5S_SCALAR);
	hid_t deltaAttributeId = H5Acreate2(subConcGroupId, "deltaTime", H5T_IEEE_F64LE,
			deltaDataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(deltaAttributeId, H5T_IEEE_F64LE, &deltaTime);
	status = H5Sclose(deltaDataspaceId);
	status = H5Aclose(deltaAttributeId);

	// Overwrite the last time step attribute of the concentration group
	hid_t lastDataspaceId = H5Screate(H5S_SCALAR);
	hid_t lastAttributeId = H5Aopen(concentrationGroupId, "lastTimeStep", H5P_DEFAULT);
	status = H5Awrite(lastAttributeId, H5T_STD_I32LE, &timeStep);
	status = H5Sclose(lastDataspaceId);
	status = H5Aclose(lastAttributeId);

	// Create property list for independent dataset write
	// (needed to be able to write the datasets without having
	// HDF5 screaming).
	propertyListId = H5Pcreate(H5P_DATASET_XFER);
	status = H5Pset_dxpl_mpio(propertyListId, H5FD_MPIO_INDEPENDENT);

	return;
}

void HDF5Utils::addConcentrationDataset(int size, int i, int j, int k) {
	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "position_" << i << "_" << j << "_" << k;

	// Create the dataspace for the dataset with dimension dims
	hsize_t dims[2];
	dims[0] = size;
	dims[1] = 2;
	concDataspaceId = H5Screate_simple(2, dims, NULL);

	// Create the dataset of concentrations for this position
	hid_t datasetId = H5Dcreate2(subConcGroupId, datasetName.str().c_str(),
	H5T_IEEE_F64LE, concDataspaceId,
	H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Create property list for independent dataset write.
	propertyListId = H5Pcreate(H5P_DATASET_XFER);
	status = H5Pset_dxpl_mpio(propertyListId, H5FD_MPIO_INDEPENDENT);

	// Close everything
	status = H5Sclose(concDataspaceId);
	status = H5Dclose(datasetId);

	return;
}

void HDF5Utils::fillConcentrations(std::vector<std::vector<double> > concVector,
		int i, int j, int k) {
	// Create the concentration array
	double concArray[concVector.size()][2];

	// Fill it with the concentration vector
	for (int i = 0; i < concVector.size(); i++) {
		concArray[i][0] = concVector.at(i).at(0);
		concArray[i][1] = concVector.at(i).at(1);
	}

	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "position_" << i << "_" << j << "_" << k;

	// Open the already created dataset of concentrations for this position
	hid_t datasetId = H5Dopen(subConcGroupId, datasetName.str().c_str(),
			H5P_DEFAULT);

	// Write concArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, propertyListId,
			&concArray);

	// Close everything
	status = H5Dclose(datasetId);

	return;
}

void HDF5Utils::finalizeFile() {
	// Close everything
	status = H5Gclose(headerGroupId);
	status = H5Gclose(networkGroupId);
	status = H5Gclose(concentrationGroupId);
	status = H5Fclose(fileId);

	return;
}

void HDF5Utils::closeFile() {
	// Close everything
	status = H5Pclose(propertyListId);
	status = H5Gclose(subConcGroupId);
	status = H5Gclose(concentrationGroupId);
	status = H5Fclose(fileId);

	return;
}

void HDF5Utils::readHeader(std::string fileName, int &nx, double &hx, int &ny,
		double &hy, int &nz, double &hz) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Open the header group
	hid_t groupId = H5Gopen(fileId, "/headerGroup", H5P_DEFAULT);

	// Open and read the nx attribute
	hid_t attributeId = H5Aopen(groupId, "nx", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_STD_I32LE, &nx);
	status = H5Aclose(attributeId);
	// Open and read the hx attribute
	attributeId = H5Aopen(groupId, "hx", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &hx);
	status = H5Aclose(attributeId);

	// Open and read the ny attribute
	attributeId = H5Aopen(groupId, "ny", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_STD_I32LE, &ny);
	status = H5Aclose(attributeId);
	// Open and read the hy attribute
	attributeId = H5Aopen(groupId, "hy", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &hy);
	status = H5Aclose(attributeId);

	// Open and read the nz attribute
	attributeId = H5Aopen(groupId, "nz", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_STD_I32LE, &nz);
	status = H5Aclose(attributeId);
	// Open and read the hz attribute
	attributeId = H5Aopen(groupId, "hz", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &hz);
	status = H5Aclose(attributeId);

	// Close everything
	status = H5Gclose(groupId);
	status = H5Fclose(fileId);

	return;
}

bool HDF5Utils::hasConcentrationGroup(std::string fileName,
		int &lastTimeStep) {
	// Initialize the boolean to return
	bool hasGroup = true;

	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Check the group
	bool groupExist = H5Lexists(fileId, "/concentrationsGroup", H5P_DEFAULT);
	// If the group exist
	if (groupExist) {
		// Open the concentration group
		concentrationGroupId = H5Gopen(fileId, "/concentrationsGroup", H5P_DEFAULT);

		// Open and read the lastTimeStep attribute
		hid_t lastAttributeId = H5Aopen(concentrationGroupId, "lastTimeStep", H5P_DEFAULT);
		status = H5Aread(lastAttributeId, H5T_STD_I32LE, &lastTimeStep);
		status = H5Aclose(lastAttributeId);

		status = H5Gclose(concentrationGroupId);

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

void HDF5Utils::readTimes(std::string fileName, int lastTimeStep, double &time,
		double &deltaTime) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentrationsGroup/concentration_" << lastTimeStep;

	// Open this specific concentration sub group
	subConcGroupId = H5Gopen(fileId, subGroupName.str().c_str(), H5P_DEFAULT);

	// Open and read the absoluteTime attribute
	hid_t timeAttributeId = H5Aopen(subConcGroupId, "absoluteTime", H5P_DEFAULT);
	status = H5Aread(timeAttributeId, H5T_IEEE_F64LE, &time);
	status = H5Aclose(timeAttributeId);

	// Open and read the deltaTime attribute
	hid_t deltaAttributeId = H5Aopen(subConcGroupId, "deltaTime", H5P_DEFAULT);
	status = H5Aread(deltaAttributeId, H5T_IEEE_F64LE, &deltaTime);
	status = H5Aclose(deltaAttributeId);

	// Close everything
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return;
}

std::vector<std::vector<double> > HDF5Utils::readNetwork(std::string fileName) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Open the dataset
	hid_t datasetId = H5Dopen(fileId, "/networkGroup/network", H5P_DEFAULT);

	// Open and read the networkSize attribute
	hid_t networkSizeAttributeId = H5Aopen(datasetId, "networkSize", H5P_DEFAULT);
	int networkSize = 0;
	status = H5Aread(networkSizeAttributeId, H5T_STD_I32LE, &networkSize);
	status = H5Aclose(networkSizeAttributeId);

	// Create the array that will receive the network
	double networkArray[networkSize][6];

	// Read the data set
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			&networkArray);

	// Fill the vector to return with the dataset
	std::vector<std::vector<double> > networkVector;
	// Loop on the size of the network
	for (int i = 0; i < networkSize; i++) {
		// Create the line to give to the vector
		std::vector<double> line;
		for (int j = 0; j < 6; j++) {
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
		int lastTimeStep, int i, int j, int k) {
	// Create te vector to return
	std::vector< std::vector<double> > toReturn;

	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "concentrationsGroup/concentration_" << lastTimeStep
			<< "/position_" << i << "_" << j << "_" << k;

	// Check the dataset
	bool datasetExist = H5Lexists(fileId, datasetName.str().c_str(),
			H5P_DEFAULT);
	// If the dataset exists
	if (datasetExist) {
		// Open the dataset
		hid_t datasetId = H5Dopen(fileId, datasetName.str().c_str(),
				H5P_DEFAULT);

		// Get the dataspace object
		hid_t dataspaceId = H5Dget_space(datasetId);

		// Get the dimensions of the dataset
		hsize_t dims[2];
		status = H5Sget_simple_extent_dims(dataspaceId, dims, H5P_DEFAULT);

		// Create the array that will receive the concentrations
		double conc[dims[0]][dims[1]];

		// Read the data set
		status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
				H5P_DEFAULT, &conc);

		// Loop on the length
		for (int i = 0; i < dims[0]; i++) {
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
