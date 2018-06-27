#include "MPIUtils.h"
#include <fstream>
#include <sstream>

// The MPI communicator
static MPI_Comm xolotlMPIComm = MPI_COMM_WORLD;

void xolotlCore::MPIUtils::initialize(MPI_Comm comm) {
	xolotlMPIComm = comm;
	
	return;
}

MPI_Comm xolotlCore::MPIUtils::getMPIComm(void) {
	return xolotlMPIComm;
}
