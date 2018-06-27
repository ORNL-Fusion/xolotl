#ifndef MPIUTILS_H
#define MPIUTILS_H

#include <mpi.h>

namespace xolotlCore {

namespace MPIUtils {

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
}

} /* namespace xolotlCore */
#endif
