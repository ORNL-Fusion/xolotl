
#include "MPIUtils.h"
#include <fstream>
#include <sstream>

using namespace xolotlCore;
using std::shared_ptr;


shared_ptr<std::istream> MPIUtils::broadcastStream(
	shared_ptr<std::istream> stream, int root, MPI_Comm comm) {
	
	int rank;
	int tasks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &tasks);
	
	shared_ptr<std::stringstream> bufferSS(new std::stringstream);
	int bufferSize;
	char *buffer;
	
	// Master task
	if (rank == root) {
		// Load the data from the input stream into memory
		
		// This method is clean and potentially bug-free, but it copies memory
		// three times. However, this is not a problem for files under
		// a few MiB, and reading from disk takes longer anyway.
		
		(*bufferSS) << stream->rdbuf();
		std::string bufferString = bufferSS->str();
		
		bufferSize = bufferString.size();
		buffer = new char[bufferSize];
		bufferString.copy(buffer, bufferSize);
	}
	
	MPI_Bcast(&bufferSize, 1, MPI_INT, root, MPI_COMM_WORLD);
	
	if (rank != root) {
		buffer = new char[bufferSize];
	}
	
	MPI_Bcast(buffer, bufferSize, MPI_CHAR, root, MPI_COMM_WORLD);
	
	// Slave tasks
	if (rank != root) {
		std::string bufferString(buffer, bufferSize);
		bufferSS->str(bufferString);
	}
	
	delete buffer;
	
	// Rewind the stringstream before returning it
	
	bufferSS->seekg(0);
	return bufferSS;
}
