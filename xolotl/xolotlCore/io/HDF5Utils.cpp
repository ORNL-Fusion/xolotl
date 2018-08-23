#include "HDF5Utils.h"
#include <PSICluster.h>
#include <iostream>
#include <sstream>
#include <hdf5.h>
#include <mpi.h>

using namespace xolotlCore;

hid_t propertyListId, fileId, concentrationGroupId, subConcGroupId,
		concDataspaceId, headerGroupId;
herr_t status;

void HDF5Utils::initializeFile(const std::string& fileName) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Create the file
	fileId = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Create the group where the header will be stored
	headerGroupId = H5Gcreate2(fileId, "headerGroup", H5P_DEFAULT, H5P_DEFAULT,
	H5P_DEFAULT);

	// Create the group where the concentrations will be stored
	concentrationGroupId = H5Gcreate2(fileId, "concentrationsGroup",
			H5P_DEFAULT,
			H5P_DEFAULT, H5P_DEFAULT);

	// Create, write, and close the last written time step attribute
	int lastTimeStep = -1;
	hid_t lastDataspaceId = H5Screate(H5S_SCALAR);
	hid_t lastAttributeId = H5Acreate2(concentrationGroupId, "lastTimeStep",
			H5T_STD_I32LE, lastDataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(lastAttributeId, H5T_STD_I32LE, &lastTimeStep);

	// Close everything
	status = H5Aclose(lastAttributeId);
	status = H5Sclose(lastDataspaceId);

	return;
}

void HDF5Utils::openFile(const std::string& fileName) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read and write access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Open the concentration group
	concentrationGroupId = H5Gopen(fileId, "concentrationsGroup", H5P_DEFAULT);
	return;
}

void HDF5Utils::fillHeader(int nx, double hx, int ny, double hy, int nz,
		double hz) {
	// Create, write, and close the nx attribute
	hid_t dataspaceId = H5Screate(H5S_SCALAR);
	hid_t attributeId = H5Acreate2(headerGroupId, "nx", H5T_STD_I32LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &nx);
	status = H5Aclose(attributeId);
	// Create, write, and close the hx attribute
	attributeId = H5Acreate2(headerGroupId, "hx", H5T_IEEE_F64LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &hx);
	status = H5Aclose(attributeId);

	// Create, write, and close the ny attribute
	attributeId = H5Acreate2(headerGroupId, "ny", H5T_STD_I32LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &ny);
	status = H5Aclose(attributeId);
	// Create, write, and close the hy attribute
	attributeId = H5Acreate2(headerGroupId, "hy", H5T_IEEE_F64LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &hy);
	status = H5Aclose(attributeId);

	// Create, write, and close the nz attribute
	attributeId = H5Acreate2(headerGroupId, "nz", H5T_STD_I32LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &nz);
	status = H5Aclose(attributeId);
	// Create, write, and close the hz attribute
	attributeId = H5Acreate2(headerGroupId, "hz", H5T_IEEE_F64LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &hz);

	// Close everything
	status = H5Aclose(attributeId);
	status = H5Sclose(dataspaceId);

	return;
}

void HDF5Utils::fillNetwork(const std::string& fileName) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	hid_t fromFileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY,
			propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	status = H5Ocopy(fromFileId, "/networkGroup", fileId, "/networkGroup",
	H5P_DEFAULT, H5P_DEFAULT);

	// Close the from file
	status = H5Fclose(fromFileId);

	return;
}

void HDF5Utils::addConcentrationSubGroup(int timeStep, double time,
		double previousTime, double deltaTime) {
	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentration_" << timeStep;

	// Create the subgroup where the concentrations at this time step will be stored
	subConcGroupId = H5Gcreate2(concentrationGroupId,
			subGroupName.str().c_str(),
			H5P_DEFAULT,
			H5P_DEFAULT, H5P_DEFAULT);

	// Create, write, and close the absolute time attribute
	hid_t dataspaceId = H5Screate(H5S_SCALAR);
	hid_t attributeId = H5Acreate2(subConcGroupId, "absoluteTime",
			H5T_IEEE_F64LE, dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &time);
	status = H5Aclose(attributeId);

	// Create, write, and close the previous time attribute
	attributeId = H5Acreate2(subConcGroupId, "previousTime", H5T_IEEE_F64LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &previousTime);
	status = H5Aclose(attributeId);

	// Create, write, and close the timestep time attribute
	attributeId = H5Acreate2(subConcGroupId, "deltaTime", H5T_IEEE_F64LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &deltaTime);
	status = H5Aclose(attributeId);
	status = H5Sclose(dataspaceId);

	// Overwrite the last time step attribute of the concentration group
	attributeId = H5Aopen(concentrationGroupId, "lastTimeStep", H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &timeStep);
	status = H5Aclose(attributeId);

	// Create property list for independent dataset write
	// (needed to be able to write the datasets without having
	// HDF5 screaming).
	propertyListId = H5Pcreate(H5P_DATASET_XFER);
	status = H5Pset_dxpl_mpio(propertyListId, H5FD_MPIO_INDEPENDENT);

	return;
}

void HDF5Utils::writeSurface1D(int timeStep, int iSurface, double nInter,
		double previousFlux) {
	// Create, write, and close the surface position attribute
	hid_t dataspaceId = H5Screate(H5S_SCALAR);
	hid_t attributeId = H5Acreate2(subConcGroupId, "iSurface", H5T_STD_I32LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &iSurface);
	status = H5Aclose(attributeId);

	// Create, write, and close the quantity of interstitial attribute
	attributeId = H5Acreate2(subConcGroupId, "nInterstitial", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &nInter);
	status = H5Aclose(attributeId);

	// Create, write, and close the flux of interstitial attribute
	attributeId = H5Acreate2(subConcGroupId, "previousIFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &previousFlux);
	status = H5Aclose(attributeId);

	// Close the dataspace
	status = H5Sclose(dataspaceId);

	return;
}

void HDF5Utils::writeSurface2D(int timeStep, std::vector<int> iSurface,
		std::vector<double> nInter, std::vector<double> previousFlux) {
	// Create the array that will store the indices and fill it
	int size = iSurface.size();
	int indexArray[size];
	for (int i = 0; i < size; i++) {
		indexArray[i] = iSurface[i];
	}

	// Create the dataspace for the dataset with dimension dims
	hsize_t dims[1];
	dims[0] = size;
	hid_t dataspaceId = H5Screate_simple(1, dims, NULL);

	// Create the dataset for the surface indices
	hid_t datasetId = H5Dcreate2(subConcGroupId, "iSurface", H5T_STD_I32LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write networkArray in the dataset
	status = H5Dwrite(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &indexArray);

	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the array that will store the quantities and fill it
	double quantityArray[size];
	for (int i = 0; i < size; i++) {
		quantityArray[i] = nInter[i];
	}

	// Create the dataset for the surface indices
	datasetId = H5Dcreate2(subConcGroupId, "nInterstitial", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write networkArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantityArray);

	// Close the dataset
	status = H5Dclose(datasetId);

	// Fill the array with the previous flux
	for (int i = 0; i < size; i++) {
		quantityArray[i] = previousFlux[i];
	}

	// Create the dataset for the surface indices
	datasetId = H5Dcreate2(subConcGroupId, "previousIFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write networkArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantityArray);

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);

	return;
}

void HDF5Utils::writeSurface3D(int timeStep,
		std::vector<std::vector<int> > iSurface,
		std::vector<std::vector<double> > nInter,
		std::vector<std::vector<double> > previousFlux) {
	// Create the array that will store the indices and fill it
	int xSize = iSurface.size();
	int ySize = iSurface[0].size();
	int indexArray[xSize][ySize];
	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {
			indexArray[i][j] = iSurface[i][j];
		}
	}

	// Create the dataspace for the dataset with dimension dims
	hsize_t dims[2];
	dims[0] = xSize;
	dims[1] = ySize;
	hid_t dataspaceId = H5Screate_simple(2, dims, NULL);

	// Create the dataset for the surface indices
	hid_t datasetId = H5Dcreate2(subConcGroupId, "iSurface", H5T_STD_I32LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write in the dataset
	status = H5Dwrite(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &indexArray);
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the array that will store the interstitial quantities and fill it
	double quantityArray[xSize][ySize];
	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {
			quantityArray[i][j] = nInter[i][j];
		}
	}

	// Create the dataset for the interstitial quantities
	datasetId = H5Dcreate2(subConcGroupId, "nInterstitial", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantityArray);
	// Close the dataset
	status = H5Dclose(datasetId);

	// Fill the array that will store the interstitial flux
	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {
			quantityArray[i][j] = previousFlux[i][j];
		}
	}

	// Create the dataset for the interstitial quantities
	datasetId = H5Dcreate2(subConcGroupId, "previousIFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantityArray);
	// Close the dataset
	status = H5Dclose(datasetId);

	// Close the dataspace
	status = H5Sclose(dataspaceId);

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

void HDF5Utils::fillConcentrations(
		const std::vector<std::vector<double> >& concVector, int i, int j,
		int k) {
	// Create the concentration array
	double concArray[concVector.size()][2];

	// Fill it with the concentration vector
	for (unsigned int n = 0; n < concVector.size(); n++) {
		concArray[n][0] = concVector.at(n).at(0);
		concArray[n][1] = concVector.at(n).at(1);
	}

	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "position_" << i << "_" << j << "_" << k;

	// Open the already created dataset of concentrations for this position
	hid_t datasetId = H5Dopen(subConcGroupId, datasetName.str().c_str(),
	H5P_DEFAULT);

	// Write concArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
			propertyListId, &concArray);

	// Close everything
	status = H5Dclose(datasetId);

	return;
}

void HDF5Utils::finalizeFile() {
	// Close everything
	status = H5Gclose(headerGroupId);
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

void HDF5Utils::readHeader(const std::string& fileName, int &nx, double &hx,
		int &ny, double &hy, int &nz, double &hz) {
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

bool HDF5Utils::hasConcentrationGroup(const std::string& fileName,
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
		concentrationGroupId = H5Gopen(fileId, "/concentrationsGroup",
				H5P_DEFAULT);

		// Open and read the lastTimeStep attribute
		hid_t lastAttributeId = H5Aopen(concentrationGroupId, "lastTimeStep",
				H5P_DEFAULT);
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

void HDF5Utils::readTimes(const std::string& fileName, int lastTimeStep,
		double &time, double &deltaTime) {
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
	hid_t attributeId = H5Aopen(subConcGroupId, "absoluteTime", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &time);
	status = H5Aclose(attributeId);

	// Open and read the deltaTime attribute
	attributeId = H5Aopen(subConcGroupId, "deltaTime", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &deltaTime);
	status = H5Aclose(attributeId);

	// Close everything
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return;
}

double HDF5Utils::readPreviousTime(const std::string& fileName,
		int lastTimeStep) {
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

	// Open and read the previousTime attribute
	double previousTime = 0.0;
	hid_t attributeId = H5Aopen(subConcGroupId, "previousTime", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &previousTime);

	// Close everything
	status = H5Aclose(attributeId);
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return previousTime;
}

int HDF5Utils::readSurface1D(const std::string& fileName, int lastTimeStep) {
	// Initialize the surface position
	int iSurface = 0;

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

	// Open and read the iSurface attribute
	hid_t attributeId = H5Aopen(subConcGroupId, "iSurface", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_STD_I32LE, &iSurface);
	status = H5Aclose(attributeId);

	// Close everything
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return iSurface;
}

std::vector<int> HDF5Utils::readSurface2D(const std::string& fileName,
		int lastTimeStep) {
	// Create the vector to return
	std::vector<int> toReturn;

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

	// Open the dataset
	hid_t datasetId = H5Dopen(subConcGroupId, "iSurface", H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
	hsize_t dims[1];
	status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);

	// Create the array that will receive the indices
	int index[dims[0]];

	// Read the data set
	status = H5Dread(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &index);

	// Loop on the length and fill the vector to return
	for (int i = 0; i < dims[0]; i++) {
		toReturn.push_back(index[i]);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return toReturn;
}

std::vector<std::vector<int> > HDF5Utils::readSurface3D(
		const std::string& fileName, int lastTimeStep) {
	// Create the vector to return
	std::vector<std::vector<int> > toReturn;

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

	// Open the dataset
	hid_t datasetId = H5Dopen(subConcGroupId, "iSurface", H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
	hsize_t dims[2];
	status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);

	// Create the array that will receive the indices
	int index[dims[0]][dims[1]];

	// Read the data set
	status = H5Dread(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &index);

	// Loop on the length and fill the vector to return
	for (int i = 0; i < dims[0]; i++) {
		// Create a temporary vector
		std::vector<int> temp;
		for (int j = 0; j < dims[1]; j++) {
			temp.push_back(index[i][j]);
		}
		// Add the temporary vector to the one to return
		toReturn.push_back(temp);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return toReturn;
}

double HDF5Utils::readNInterstitial1D(const std::string& fileName,
		int lastTimeStep) {
	// Initialize the surface position
	double nInter = 0.0;

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

	// Open and read the iSurface attribute
	hid_t attributeId = H5Aopen(subConcGroupId, "nInterstitial", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &nInter);
	status = H5Aclose(attributeId);

	// Close everything
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return nInter;
}

std::vector<double> HDF5Utils::readNInterstitial2D(const std::string& fileName,
		int lastTimeStep) {
	// Create the vector to return
	std::vector<double> toReturn;

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

	// Open the dataset
	hid_t datasetId = H5Dopen(subConcGroupId, "nInterstitial", H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
	hsize_t dims[1];
	status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);

	// Create the array that will receive the indices
	double quantity[dims[0]];

	// Read the data set
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantity);

	// Loop on the length and fill the vector to return
	for (int i = 0; i < dims[0]; i++) {
		toReturn.push_back(quantity[i]);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return toReturn;
}

std::vector<std::vector<double> > HDF5Utils::readNInterstitial3D(
		const std::string& fileName, int lastTimeStep) {
	// Create the vector to return
	std::vector<std::vector<double> > toReturn;

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

	// Open the dataset
	hid_t datasetId = H5Dopen(subConcGroupId, "nInterstitial", H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
	hsize_t dims[2];
	status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);

	// Create the array that will receive the indices
	double quantity[dims[0]][dims[1]];

	// Read the data set
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantity);

	// Loop on the length and fill the vector to return
	for (int i = 0; i < dims[0]; i++) {
		// Create a temporary vector
		std::vector<double> temp;
		for (int j = 0; j < dims[1]; j++) {
			temp.push_back(quantity[i][j]);
		}
		// Add the temporary vector to the one to return
		toReturn.push_back(temp);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return toReturn;
}

double HDF5Utils::readPreviousIFlux1D(const std::string& fileName,
		int lastTimeStep) {
	// Initialize the surface position
	double previousFlux = 0.0;

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

	// Open and read the iSurface attribute
	hid_t attributeId = H5Aopen(subConcGroupId, "previousIFlux", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &previousFlux);
	status = H5Aclose(attributeId);

	// Close everything
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return previousFlux;
}

std::vector<double> HDF5Utils::readPreviousIFlux2D(const std::string& fileName,
		int lastTimeStep) {
	// Create the vector to return
	std::vector<double> toReturn;

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

	// Open the dataset
	hid_t datasetId = H5Dopen(subConcGroupId, "previousIFlux", H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
	hsize_t dims[1];
	status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);

	// Create the array that will receive the indices
	double flux[dims[0]];

	// Read the data set
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &flux);

	// Loop on the length and fill the vector to return
	for (int i = 0; i < dims[0]; i++) {
		toReturn.push_back(flux[i]);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return toReturn;
}

std::vector<std::vector<double> > HDF5Utils::readPreviousIFlux3D(
		const std::string& fileName, int lastTimeStep) {
	// Create the vector to return
	std::vector<std::vector<double> > toReturn;

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

	// Open the dataset
	hid_t datasetId = H5Dopen(subConcGroupId, "previousIFlux", H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
	hsize_t dims[2];
	status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);

	// Create the array that will receive the indices
	double quantity[dims[0]][dims[1]];

	// Read the data set
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantity);

	// Loop on the length and fill the vector to return
	for (int i = 0; i < dims[0]; i++) {
		// Create a temporary vector
		std::vector<double> temp;
		for (int j = 0; j < dims[1]; j++) {
			temp.push_back(quantity[i][j]);
		}
		// Add the temporary vector to the one to return
		toReturn.push_back(temp);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return toReturn;
}

std::vector<std::vector<double> > HDF5Utils::readNetwork(
		const std::string& fileName) {
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
	hid_t networkSizeAttributeId = H5Aopen(datasetId, "networkSize",
			H5P_DEFAULT);
	int networkSize = 0;
	status = H5Aread(networkSizeAttributeId, H5T_STD_I32LE, &networkSize);
	status = H5Aclose(networkSizeAttributeId);

	// Create the array that will receive the network
	double *networkArray = new double[networkSize * 6];

	// Read the data set
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			networkArray);

	// Fill the vector to return with the dataset
	std::vector<std::vector<double> > networkVector;
	// Loop on the size of the network
	for (int i = 0; i < networkSize; i++) {
		// Create the line to give to the vector
		std::vector<double> line;
		for (int j = 0; j < 6; j++) {
			line.push_back(networkArray[i * 6 + j]);
		}
		networkVector.push_back(line);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Fclose(fileId);
	delete[] networkArray;

	return networkVector;
}

std::vector<std::vector<double> > HDF5Utils::readGridPoint(
		const std::string& fileName, int lastTimeStep, int i, int j, int k) {
	// Create the vector to return
	std::vector<std::vector<double> > toReturn;

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
		status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);

		// Create the array that will receive the concentrations
		double conc[dims[0]][dims[1]];

		// Read the data set
		status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
		H5P_DEFAULT, &conc);

		// Loop on the length
		for (unsigned int n = 0; n < dims[0]; n++) {
			// Create the concentration vector for this cluster
			std::vector<double> tmp;
			tmp.push_back(conc[n][0]);
			tmp.push_back(conc[n][1]);

			// Add it to the main vector
			toReturn.push_back(tmp);
		}

		// Close everything
		status = H5Dclose(datasetId);
		status = H5Sclose(dataspaceId);
	}

	// Close the file
	status = H5Fclose(fileId);

	return toReturn;
}
