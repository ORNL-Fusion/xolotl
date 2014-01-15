
#include "MPIUtils.h"
#include <fstream>
#include <sstream>

using namespace xolotlCore;
using std::shared_ptr;


shared_ptr<std::istream> MPIUtils::broadcastStream(
	shared_ptr<std::istream> stream, int master, MPI_Comm comm) {
	
	int rank;
	int tasks;
	auto bufferSS = std::make_shared<std::stringstream>();
	int bufferSize;
	char *buffer;

	// Get the rank and number of tasks
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &tasks);
	
	// Master task
	if (rank == master) {
		// Load the data from the input stream into memory
		// This method is clean and potentially bug-free, but it copies memory
		// three times. However, this is not a problem for files under
		// a few MiB, and reading from disk takes longer anyway.
		(*bufferSS) << stream->rdbuf();
		std::string bufferString = bufferSS->str();
		// Create the new buffer and copy its contents into the string
		bufferSize = bufferString.size();
		buffer = new char[bufferSize];
		bufferString.copy(buffer, bufferSize);
	}
	
	// Broadcast the size
	MPI_Bcast(&bufferSize, 1, MPI_INT, master, MPI_COMM_WORLD);
	
	// Create a new buffer on the nodes
	if (rank != master) {
		buffer = new char[bufferSize];
	}
	
	// Broadcast the buffer
	MPI_Bcast(buffer, bufferSize, MPI_CHAR, master, MPI_COMM_WORLD);
	
	// Creat the string from the buffer on the nodes
	if (rank != master) {
		std::string bufferString(buffer, bufferSize);
		bufferSS->str(bufferString);
	}
	
	// Rewind the stringstream before returning it
	bufferSS->seekg(0);

	// Clean up the buffer memory and return it
	delete buffer;
	
	return bufferSS;
}
