#ifndef MPIBROADCASTER_H
#define MPIBROADCASTER_H

#include <mpi.h>
#include <memory>
#include <iostream>


namespace xolotlCore {

namespace MPIUtils {

	/**
	 * Sends the input buffer from the master task to all the slaves
	 *
	 * This method blocks until it is called by all processes.
	 * The root process must have a valid network stream, but the streams
	 * of worker tasks are ignored.
	 *
	 * @param stream The stream to broadcast
	 * @param root The rank of the master process
	 * @return The buffer
	 */
	std::shared_ptr<std::istream> broadcastStream(
		std::shared_ptr<std::istream> stream, int root);
}

} /* namespace xolotlCore */
#endif
