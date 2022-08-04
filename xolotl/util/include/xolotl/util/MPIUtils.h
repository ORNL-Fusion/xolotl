#ifndef MPIBROADCASTER_H
#define MPIBROADCASTER_H

#include <iostream>
#include <memory>

#include <mpi.h>

namespace xolotl
{
namespace util
{
/**
 * Specify the MPI communicator
 *
 * @param comm The communicator we want to use.
 */
void
setMPIComm(MPI_Comm comm);

/**
 * Access the MPI communicator.
 *
 *  @return The communicator
 */
MPI_Comm
getMPIComm();

/**
 * Initialize MPI with const char array
 */
void
mpiInit(int& argc, const char* argv[]);

/**
 * Has MPI been initialized
 */
bool
mpiInitialized();

/**
 * Get rank of current process within Xolotl's communicator
 */
int
getMPIRank();
} // namespace util
} // namespace xolotl
#endif
