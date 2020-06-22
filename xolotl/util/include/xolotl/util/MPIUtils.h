#ifndef MPIBROADCASTER_H
#define MPIBROADCASTER_H

#include <mpi.h>
#include <memory>
#include <iostream>


namespace xolotl
{
namespace util
{

/**
 * Initialize the MPI communicator
 *
 * @param comm The communicator we want to use.
 */
void initialize(MPI_Comm comm);

/**
 * Access the MPI communicator.
 *
 *  @return The communicator
 */
MPI_Comm getMPIComm(void);

} // namespace util
} // namespace xolotl
#endif
