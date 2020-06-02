#ifndef XCORE_DOINORDER_H
#define XCORE_DOINORDER_H

namespace xolotl {
namespace core {


/**
 * Have each MPI rank execute the given function in rank order.
 * Can pass a zero-argument functor or lambda as the function.
 * Intended to support printf-style debugging when a true parallel 
 * debugger is not an option.
 *
 * @param func The function each rank should execute.
 * @param comm The MPI communicator that defines the ranks and rank ordering.
 * @param msg A message that rank 0 should write before executing the function.
 *              No message is written if msg is empty.
 * @param os The output stream to use for writing the message, if a message 
 *              was given.
 */
template<typename F>
void
DoInOrder(F func,
            MPI_Comm comm = MPI_COMM_WORLD,
            std::string msg = "",
            std::ostream& os = std::cout) {

    int token = 7;  // value doesn't really matter
    int tag = 9;    // value doesn't really matter

    int rank = -1;
    int size = -1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPI_Barrier(comm);
    if(rank == 0) {
        if(not msg.empty()) {
            os << msg << std::endl;
        }
        func();
        MPI_Send(&token, 1, MPI_INT, (rank + 1) % size, tag, comm);
        MPI_Recv(&token, 1, MPI_INT, (size - 1), tag, comm, MPI_STATUS_IGNORE);
    }
    else {
        MPI_Recv(&token, 1, MPI_INT, (rank - 1), tag, comm, MPI_STATUS_IGNORE);
        func();
        MPI_Send(&token, 1, MPI_INT, (rank + 1) % size, tag, comm);
    }
}

} /* end namespace core */
} /* end namespace xolotl */

#endif // XCORE_DOINORDER_H
