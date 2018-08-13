#ifndef XSOLVER_MONITOR_H
#define XSOLVER_MONITOR_H

// Includes
#include <IReactionNetwork.h>

namespace xolotlSolver {

/**
 * Copy the network group (if it exists) from one file to another,
 * or write it from scratch.
 * For performance reasons, does this using a single-process communicator.
 * Files must be closed before calling this function.
 *
 * @param _comm The MPI communicator to use to determine which process
 *              should do the copy.
 * @param srcFileName The path to the file containing the network group.
 * @param targetFileName The path to the file into which the network group
 *                  should be copied.
 * @param network The network to write.
 */
void writeNetwork(MPI_Comm _comm, std::string srcFileName,
		std::string targetFileName, IReactionNetwork& network);

} // namespace xolotlSolver

#endif // XSOLVER_MONITOR_H
