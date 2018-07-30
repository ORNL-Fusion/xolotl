#ifndef XCORE_HDF5FILE_DATASET_H
#define XCORE_HDF5FILE_DATASET_H

#include <numeric>
#include "boost/range/counting_range.hpp"

namespace xolotlCore
{

template<typename T>
HDF5File::DataSetTBase<T>::DataSetTBase(const HDF5Object& loc,
                                    std::string dsetName,
                                    const DataSpace& dspace)
  : DataSetBase(loc, dsetName)
{
    setId(H5Dcreate(loc.getId(),
                        dsetName.c_str(),
                        TypeInFile<T>().getId(),
                        dspace.getId(),
                        H5P_DEFAULT,
                        H5P_DEFAULT,
                        H5P_DEFAULT));
    if(getId() < 0)
    {
        std::ostringstream estr;
        estr << "Failed to create dataset " << getName();
        throw HDF5Exception(estr.str());
    }
}

template<typename T>
HDF5File::DataSetTBase<T>::DataSetTBase(const HDF5Object& loc,
                                        std::string dsetName)
  : DataSetBase(loc, dsetName)
{
    setId(H5Dopen(loc.getId(), dsetName.c_str(), H5P_DEFAULT));
    if(getId() < 0)
    {
        std::ostringstream estr;
        estr << "Failed to open dataset " << getName();
        throw HDF5Exception(estr.str());
    }
}


// Partial specialization for vector of scalars.  Should not be used for
// vector of vectors - need a different specialization for those.
template<typename T>
std::vector<T>
HDF5File::DataSet<std::vector<T>>::read(void) const {

    // Ensure we have space for the data.
    SimpleDataSpace<1> dspace(*this);
    auto nItems = dspace.getDims()[0];
    std::vector<T> ret(nItems);

    // Read data from the file.
    TypeInMemory<T> memType;
    auto status = H5Dread(this->getId(), memType.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, ret.data());
    if(status < 0) {
        std::ostringstream estr;
        estr << "Unable to read value of vector dataset " << this->getName();
        throw HDF5Exception(estr.str());
    }
    return ret;
}


#if READY
template<typename T>
void
HDF5File::DataSet<std::vector<T>>::write(const std::vector<T>& data) const {

    // Verify the dimensionality of the attribute.
    SimpleDataSpace<1> dspace(*this);
    auto nItems = dspace.getDims()[0];
    if(nItems != data.size())
    {
        throw HDF5Exception("Size mismatch when setting vector-valued dataset");
    }

    // Write the data.
    TypeInMemory<T> memType;
	auto status = H5Dwrite(id,
                            memType.getId(),
                            H5S_ALL,
                            H5S_ALL,
	                        H5P_DEFAULT,
                            data.data());
    if(status < 0)
    {
        std::ostringstream estr;
        estr << "Failed to write variable length vector data to dataset " << name;
        throw HDF5Exception(estr.str());
    }
}


template<typename T>
void
HDF5File::VectorsDataSet<T>::write(const std::vector<std::vector<T>>& data) const
{
    // Define variable-length data we will write.
    auto nDims = data.size();
    std::vector<hvl_t> vlData(nDims);
    for(auto i : boost::counting_range<uint32_t>(0, nDims))
    {
        vlData[i].p = const_cast<void*>(static_cast<const void*>(data[i].data()));
        vlData[i].len = data[i].size();
    }

    // Write the data to the file.
    TypeInMemory<std::vector<T>> memType;
    auto status = H5Dwrite(DataSet<std::vector<T>>::id,
                        memType.getId(),
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        vlData.data());
    if(status < 0)
    {
        throw HDF5Exception("Failed to write variable length vector data to dataset");
    }

    // We do *not* call H5vlen_reclaim here, 
    // because the data buffers were owned by the vectors passed
    // in by the caller.
}

template<typename T>
std::vector<std::vector<T>>
HDF5File::VectorsDataSet<T>::read(void) const
{
    // Determine the dataspace for our data set.
    // TODO is this really a dataspace with rank 1?
    SimpleDataSpace<1> dspace(*this);
    auto nVectors = dspace.getDims()[0];

    // Read the data from the file.
    std::vector<hvl_t> vlData(nVectors);
    TypeInMemory<std::vector<T>> memType;
    auto status = H5Dread(DataSet<std::vector<T>>::id,
                            memType.getId(),
                            H5S_ALL,
                            H5S_ALL,
                            H5P_DEFAULT,
                            vlData.data());
    if(status < 0)
    {
        std::ostringstream estr;
        estr << "Failed to read vector of vector data from dataset " 
            << this->getName();
        throw HDF5Exception(estr.str());
    }

    // Convert to our output format.
    std::vector<std::vector<T>> ret(nVectors);
    for(auto i : boost::counting_range<uint32_t>(0, nVectors))
    {
        for(auto j : boost::counting_range<uint32_t>(0, vlData[i].len))
        {
            T* currData = static_cast<T*>(vlData[i].p);
            ret[i].emplace_back(currData[j]);
        }
    }

    // Release memory allocated by the HDF5 library.
    H5Dvlen_reclaim(memType.getId(),
                    dspace.getId(),
                    H5P_DEFAULT,
                    vlData.data());

    return ret;
}
#endif // READY



template<typename T>
HDF5File::RaggedDataSet2D<T>::RaggedDataSet2D(MPI_Comm _comm,
                                    const HDF5Object& loc,
                                    std::string dsetName,
                                    int baseX,
                                    const Ragged2DType& data)
  : RaggedDataSetBase(_comm),
    DataSetTBase<T>(loc, dsetName, *(buildDataSpace(_comm, data))) {

    // We assume the gridpoint values are indices into the gridpoint array,
    // so non-negative and base 0.
    assert(baseX >= 0);
    
    // Write the indexing metadata describing our part of the flattened dataset.
    uint32_t globalBaseIdx;
    uint32_t myNumItems;
    std::tie(globalBaseIdx, myNumItems) = writeStartingIndices(baseX, data);

    // Finally, write our part of the data itself.
    writeData(globalBaseIdx, myNumItems, data);
}


template<typename T>
HDF5File::RaggedDataSet2D<T>::RaggedDataSet2D(MPI_Comm _comm,
                                    const HDF5Object& loc,
                                    std::string dsetName)
  : RaggedDataSetBase(_comm),
    DataSetTBase<T>(loc, dsetName) {

    // Nothing else to do.
}



template<typename T>
std::vector<uint32_t>
HDF5File::RaggedDataSet2D<T>::findNumItemsByPoint(const Ragged2DType& data) {

    std::vector<uint32_t> ret(data.size());
    for(auto i = 0; i < data.size(); ++i) {
        ret[i] = data[i].size();
    }
    return ret;
}



template<typename T>
std::unique_ptr<HDF5File::SimpleDataSpace<1>>
HDF5File::RaggedDataSet2D<T>::buildDataSpace(MPI_Comm _comm,
                                        const Ragged2DType& data) {

    // Build the data space for the flattened data itself.
    // When a file is opened for parallel access, HDF5 seems to require 
    // that all processes make structural changes like defining a 
    // dataset as a collective operation.
    //
    // The dataspace of the flattened data is a 1D array with size
    // equal to the total number of data items across all processes.
    // The given 'data' variable only contains the data items for 
    // our own process, so we need to aggregate those values across 
    // all processes.

    // Figure out the number of items associated with each grid point.
    // TODO is there much benefit to caching this info since we 
    // also use it when figuring out starting indices?
    auto myNumItemsByPoint = findNumItemsByPoint(data);
    auto myNumPoints = data.size();
    assert(myNumItemsByPoint.size() == myNumPoints);

    // Determine the total number of items we own.
    uint32_t myNumItems = 
        std::accumulate(myNumItemsByPoint.begin(), myNumItemsByPoint.end(), 0);

    // Determine the total number of items across all processes.
    // This is the size of the 1D flattened data set.
    uint32_t totalNumItems;
    MPI_Allreduce(&myNumItems,
                    &totalNumItems,
                    1,
                    MPI_UNSIGNED,
                    MPI_SUM,
                    _comm);

    // Define the dataspace for the flattened data.
    SimpleDataSpace<1>::Dimensions dims { totalNumItems };
    std::unique_ptr<SimpleDataSpace<1>> dataspace(new SimpleDataSpace<1>(dims));
    return std::move(dataspace);
}

template<typename F>
void
DoInOrder(std::string msg, int rank, int size, F func) {
    int token = 7;
    int tag = 9;
    if(rank == 0) {
        std::cout << msg << std::endl;
        func();
        MPI_Send(&token, 1, MPI_INT, (rank + 1) % size, tag, MPI_COMM_WORLD);
        MPI_Recv(&token, 1, MPI_INT, (size - 1), tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        MPI_Recv(&token, 1, MPI_INT, (rank - 1), tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        func();
        MPI_Send(&token, 1, MPI_INT, (rank + 1) % size, tag, MPI_COMM_WORLD);
    }
}

inline
std::ostream&
operator<<(std::ostream& os, const std::pair<int, double>& p) {

    os << '(' << p.first << ',' << p.second << ')';
    return os;
}


template<typename T>
std::pair<uint32_t, uint32_t>
HDF5File::RaggedDataSet2D<T>::writeStartingIndices(int baseX, 
                                            const Ragged2DType& data) const {

    // Determine our position within the MPI communicator used to
    // access the file.
    int commRank;
    MPI_Comm_rank(comm, &commRank);
    int commSize;
    MPI_Comm_size(comm, &commSize);

#if READY
    DoInOrder("Data", commRank, commSize,
        [baseX, &data]() {
            for(auto i = 0; i < data.size(); ++i) {
                std::cout << baseX + i << ": ";
                for(auto const& p : data[i]) {
                    std::cout << ' ' << p;
                }
                std::cout << '\n';
            }
        });
#endif // READY

    // Determine the local starting indices of the data we own.
    // First, determine the number of items per grid point...
    auto myNumPoints = data.size();
    auto myNumItemsByPoint = findNumItemsByPoint(data);
    assert(myNumItemsByPoint.size() == myNumPoints);

#if READY
    MPI_Barrier(MPI_COMM_WORLD);
    DoInOrder("NumItemsByPoint", commRank, commSize,
        [commRank, &myNumItemsByPoint]() {
            std::cout << commRank << ": ";
            for(auto const& p : myNumItemsByPoint) {
                std::cout << ' ' << p;
            }
            std::cout << '\n';
        });
#endif // READY

    // ...then determine starting indices within our local flat array.
    // Note: we're using C++11's std::partial_sum which does an inclusive
    // scan, so we need to manually shift the values in the output
    // vector, and then drop the last value after the scan.
    // If we had C++17 support, we could just use std::exclusive_scan.
    FlatStartingIndicesType myStartingIndices(myNumPoints+1);
    myStartingIndices[0] = 0;
    std::partial_sum(myNumItemsByPoint.begin(), myNumItemsByPoint.end(),
                        myStartingIndices.begin() + 1);
    myStartingIndices.resize(myNumPoints);

#if READY
    MPI_Barrier(MPI_COMM_WORLD);
    DoInOrder("Local starting indices", commRank, commSize,
        [commRank, &myStartingIndices]() {
            std::cout << commRank << ": ";
            for(auto const& p : myStartingIndices) {
                std::cout << ' ' << p;
            }
            std::cout << '\n';
        });
#endif // READY

    // Convert local starting indices into global starting indices
    // for the gridpoints we own.
    uint32_t myNumItems = myStartingIndices.back() + myNumItemsByPoint.back();
    uint32_t globalBaseIdx = 0;
    MPI_Exscan(&myNumItems,
                &globalBaseIdx,
                1,
                MPI_UNSIGNED,
                MPI_SUM,
                comm);
    if(commRank == 0) {
        // Ensure global base index is 0 in rank 0.
        // Some (all?) MPI_Exscan implementations leave rank 0 output
        // value undefined.
        globalBaseIdx = 0;
    }
    FlatStartingIndicesType globalStartingIndices(myStartingIndices.size());
    std::transform(myStartingIndices.begin(), myStartingIndices.end(),
                    globalStartingIndices.begin(),
                    [globalBaseIdx](uint32_t idx) -> uint32_t {
                        return idx + globalBaseIdx;
                    });
    assert(globalStartingIndices.size() == myNumPoints);

    //
    // Write the starting index metadata.
    // Instead of keeping a separate dataset for the number of
    // items for each grid point, we only keep the starting index
    // for each grid point's items.  Then we determine the number of 
    // items for grid point i as (startingIndices[i+1] - startingIndices[i]).
    // To avoid having a special case for the last grid point, we
    // extend the startingIndices dataset with one element indicating
    // the total number of items.
    //

    // Determine the total number of points across all processes.
    uint32_t totalNumPoints;
    MPI_Allreduce(&myNumPoints,
                    &totalNumPoints,
                    1,
                    MPI_UNSIGNED,
                    MPI_SUM,
                    comm);

    // Create the global index dataspace.
    // Recall that we make it one larger than the total number of
    // grid points so that we can store the total number of items 
    // and use the (i+1) - i approach.
    SimpleDataSpace<1>::Dimensions globalIndexDims { totalNumPoints + 1 };
    SimpleDataSpace<1> indexDataSpace(globalIndexDims);

    // Create the global index dataset.
    std::ostringstream indexDatasetNameStr;
    indexDatasetNameStr << getName() << startIndicesDatasetNameSuffix;
    DataSet<uint32_t> indexDataset(this->getLocation(),
                                    indexDatasetNameStr.str(),
                                    indexDataSpace);

    // Write the index dataset in parallel.
    //
    // Describe our part of the data within the global indexing dataspace.
    SimpleDataSpace<1>::Dimensions indexCounts { myNumPoints };
    SimpleDataSpace<1>::Dimensions indexOffsets { (uint32_t)baseX };
    if(commRank == (commSize - 1)) {
        // We are the last MPI rank.  Add the total number of items
        // to our part of the globalStartingIndices dataset.
        ++indexCounts[0];
        SimpleDataSpace<1> globalDataspace(*this);
        auto globalDims = globalDataspace.getDims();
        auto totalNumItems = globalDims[0];
        globalStartingIndices.push_back(totalNumItems);
    }
    SimpleDataSpace<1> indexMemspace(indexCounts);
#if READY
    MPI_Barrier(MPI_COMM_WORLD);
    DoInOrder("Global starting indices", commRank, commSize,
        [commRank, &globalStartingIndices]() {
            std::cout << commRank << ": ";
            for(auto const& p : globalStartingIndices) {
                std::cout << ' ' << p;
            }
            std::cout << '\n';
        });
#endif // READY

#if READY
    MPI_Barrier(MPI_COMM_WORLD);
    DoInOrder("Index array parts", commRank, commSize,
        [commRank, baseX, &indexCounts]() {
            std::cout << commRank << ": "
                << '[' << baseX << ", " << (baseX + indexCounts[0]) << ')'
                << std::endl;
        });
#endif // READY


    // Select our hyperslab within the file.
    SimpleDataSpace<1> indexFilespace(indexDataset);
    H5Sselect_hyperslab(indexFilespace.getId(),
                        H5S_SELECT_SET,
                        indexOffsets.data(),
                        nullptr,
                        indexCounts.data(),
                        nullptr);

    // Write the index metadata using a collective write.
    PropertyList plist(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist.getId(), H5FD_MPIO_COLLECTIVE);
    TypeInMemory<uint32_t> indexMemType;
    H5Dwrite(indexDataset.getId(),
                indexMemType.getId(),
                indexMemspace.getId(),
                indexFilespace.getId(),
                plist.getId(),
                globalStartingIndices.data());

    return std::make_pair(globalBaseIdx, myNumItems);
}

template<typename T>
void
HDF5File::RaggedDataSet2D<T>::writeData(uint32_t globalBaseIdx, 
                                        uint32_t myNumItems,
                                        const Ragged2DType& data) const {

    // Describe our data within the global dataspace.
    SimpleDataSpace<1>::Dimensions dataCounts {myNumItems};
    SimpleDataSpace<1>::Dimensions dataOffsets {globalBaseIdx};
    SimpleDataSpace<1> dataMemSpace(dataCounts);

    // Select our hyperslab within the file.
    SimpleDataSpace<1> dataFileSpace(*this);
    auto status = H5Sselect_hyperslab(dataFileSpace.getId(),
                                        H5S_SELECT_SET,
                                        dataOffsets.data(),
                                        nullptr,
                                        dataCounts.data(),
                                        nullptr);
    if(status < 0) {
        std::ostringstream estr;
        estr << "Failed to select our part of dataset " << getName();
        throw HDF5Exception(estr.str());
    }

    // Flatten our data.
    FlatType flatData;
    for(auto const& currItems : data) {
        std::copy(currItems.begin(), currItems.end(),
                std::back_inserter(flatData));
    }
    if(flatData.size() != myNumItems) {
        std::cerr << "flatData.size() = " << flatData.size()
            << ", myNumItems = " << myNumItems
            << std::endl;
    }
    assert(flatData.size() == myNumItems);

    // Write the flat data using a collective write.
    PropertyList plist(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist.getId(), H5FD_MPIO_COLLECTIVE);
    TypeInMemory<T> memType;
    status = H5Dwrite(getId(),
                memType.getId(),
                dataMemSpace.getId(),
                dataFileSpace.getId(),
                plist.getId(),
                flatData.data());
    if(status < 0)
    {
        std::ostringstream estr;
        estr << "Failed to write dataset " << getName();
        throw HDF5Exception(estr.str());
    }
}


template<typename T>
typename HDF5File::RaggedDataSet2D<T>::Ragged2DType
HDF5File::RaggedDataSet2D<T>::read(int baseX, int numX) const {

    // We assume the gridpoint values are indices into the gridpoint array,
    // so non-negative and base 0.
    assert(baseX >= 0);
    
    // Read our part of the indexing metadata.
    auto globalStartingIndices = readStartingIndices(baseX, numX);

    // Read our part of the actual data.
    auto ret = readData(globalStartingIndices);
#if READY
    int commRank;
    int commSize;
    MPI_Comm_rank(comm, &commRank);
    MPI_Comm_size(comm, &commSize);

    MPI_Barrier(comm);
    DoInOrder("Read ragged data", commRank, commSize,
            [baseX, numX, &ret]() {
                for(auto i = 0; i < numX; ++i ) {
                    std::cout << baseX + i << ':';
                    for(const auto& currData : ret[i]) {
                        std::cout << ' ' << currData;
                    }
                    std::cout << std::endl;
                }
            });
#endif // READY

    return ret;
}



template<typename T>
std::vector<uint32_t>
HDF5File::RaggedDataSet2D<T>::readStartingIndices(int baseX, int numX) const {

    // Determine our position within the MPI communicator used to
    // access the file.
    int commRank;
    MPI_Comm_rank(comm, &commRank);
    int commSize;
    MPI_Comm_size(comm, &commSize);

    // Open our index dataset.
    std::ostringstream indexDatasetNameStr;
    indexDatasetNameStr << getName() << startIndicesDatasetNameSuffix;
    DataSet<uint32_t> indexDataset(this->getLocation(),
                                    indexDatasetNameStr.str());

    // Determine how many index values we are going to read.
    // and how many of those that we own.
    // Unlike during a write where only one rank extends the number
    // of index values it writes, when we read *each* rank reads
    // one extra so that it can easily compute how many data items
    // it owns.
    auto myNumIndexValues = numX + 1;
    SimpleDataSpace<1>::Dimensions indexCounts { (uint32_t)myNumIndexValues };
    SimpleDataSpace<1>::Dimensions indexOffsets { (uint32_t)baseX };
    SimpleDataSpace<1> indexMemspace(indexCounts);


    // Select our hyperslab within the file.
    SimpleDataSpace<1> indexFilespace(indexDataset);
    H5Sselect_hyperslab(indexFilespace.getId(),
                        H5S_SELECT_SET,
                        indexOffsets.data(),
                        nullptr,
                        indexCounts.data(),
                        nullptr);

    // Read the index metadata we need using a collective operation.
    std::vector<uint32_t> globalStartingIndices(myNumIndexValues);
    PropertyList plist(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist.getId(), H5FD_MPIO_COLLECTIVE);
    TypeInMemory<uint32_t> indexMemType;
    H5Dread(indexDataset.getId(),
                indexMemType.getId(),
                indexMemspace.getId(),
                indexFilespace.getId(),
                plist.getId(),
                globalStartingIndices.data());

#if READY
    MPI_Barrier(comm);
    DoInOrder("global starting indices", commRank, commSize,
            [commRank, &globalStartingIndices]() {
                std::cout << commRank << ": ";
                for(auto const& p : globalStartingIndices) {
                    std::cout << ' ' << p;
                }
                std::cout << '\n';
            });
#endif // READY

    return globalStartingIndices;
}


template<typename T>
typename HDF5File::RaggedDataSet2D<T>::Ragged2DType
HDF5File::RaggedDataSet2D<T>::readData(
                    const std::vector<uint32_t>& globalStartingIndices) const {

#if READY
    // Determine our position within the MPI communicator used to
    // access the file.
    int commRank;
    MPI_Comm_rank(comm, &commRank);
    int commSize;
    MPI_Comm_size(comm, &commSize);
#endif // READY

    // Determine the total number of values we own.
    auto myNumItems = globalStartingIndices.back() - globalStartingIndices[0];
    
    // Specify our part of the dataset to read.
    SimpleDataSpace<1>::Dimensions dataCounts { (uint32_t)myNumItems };
    SimpleDataSpace<1>::Dimensions dataOffsets { (uint32_t)globalStartingIndices[0] };
    SimpleDataSpace<1> dataMemspace(dataCounts);

    // Select our hyperslab within the file.
    SimpleDataSpace<1> dataFilespace(*this);
    H5Sselect_hyperslab(dataFilespace.getId(),
                        H5S_SELECT_SET,
                        dataOffsets.data(),
                        nullptr,
                        dataCounts.data(),
                        nullptr);

    // Read the index metadata we need using a collective operation.
    FlatType flatData(myNumItems);
    PropertyList plist(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist.getId(), H5FD_MPIO_COLLECTIVE);
    TypeInMemory<T> dataMemType;
    H5Dread(getId(),
                dataMemType.getId(),
                dataMemspace.getId(),
                dataFilespace.getId(),
                plist.getId(),
                flatData.data());
#if READY
    MPI_Barrier(comm);
    DoInOrder("flatData", commRank, commSize,
            [commRank, &flatData]() {
                std::cout << commRank << ": ";
                for(auto const& d : flatData){
                    std::cout << ' ' << d;
                }
                std::cout << '\n';
            });
#endif // READY

    // Convert from flat representation to ragged 2D representation.
    //
    // First, determine number of points we own.
    // Our number of points owned is same as number of starting indices,
    // except we reduce by 1 to remove the one extra that we read
    // to help us compute the number of values we own.
    auto myNumPoints = globalStartingIndices.size() - 1;

    // Next, compute local starting indices from the global starting
    // indices we have.
    std::vector<uint32_t> localStartingIndices(globalStartingIndices.size());
    uint32_t globalBaseIdx = globalStartingIndices[0];
    std::transform(globalStartingIndices.begin(),
                    globalStartingIndices.end(),
                    localStartingIndices.begin(),
                    [globalBaseIdx](uint32_t v) -> uint32_t {
                        return v - globalBaseIdx;
                    });

    Ragged2DType ret(myNumPoints);
    for(auto ptIdx = 0; ptIdx < myNumPoints; ++ptIdx) {
        auto localBaseIdx = localStartingIndices[ptIdx];
        auto numValues = localStartingIndices[ptIdx + 1] - localBaseIdx;
        ret[ptIdx].reserve(numValues);
        for(auto valIdx = 0; valIdx < numValues; ++valIdx) {
            ret[ptIdx].emplace_back(flatData[localBaseIdx + valIdx]);
        }
    }

    return ret;
}


} // namespace xolotlCore

#endif // XCORE_HDF5FILE_DATASET_H
