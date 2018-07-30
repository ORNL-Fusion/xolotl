#ifndef XCORE_HDF5FILE_GROUP_H
#define XCORE_HDF5FILE_GROUP_H

namespace xolotlCore {

#if READY
template<typename T>
void HDF5File::Group::write2DRaggedDataset(MPI_Comm comm,
                                    int baseX,
                                    const Ragged2DType<T>& raggedData) const {

    // Determine our place in the MPI communicator used for accessing our file.
    int commSize = -1;
    MPI_Comm_size(comm, &commSize);
    int commRank = -1;
    MPI_Comm_rank(comm, &commRank);

    // Write the indexing dataset and the flattened data.
    writeStartingIndexDataset<T>(comm, commRank, commSize, baseX, raggedData);
#if READY
    writeFlattenedDataset<T>(comm, commRank, commSize, baseX, raggedData);
#endif // READY
}


template<typename T>
void HDF5File::Group::writeStartingIndexDataset(MPI_Comm comm,
                            int commRank,
                            int commSize,
                            int baseX,
                            const Ragged2DType<T>& data) const {


    // Determine the starting index of the concentrations for
    // each of our gridpoints within the global, flattened data space.
    //
    // First, figure out how many items we have for each x value...
    auto nLocalPoints = data.size();
    FlatStartingIndicesType nItemsByPoint;
    for(auto i = 0; i < nLocalPoints; ++i) {
        nItemsByPoint[i] = data[i].size();
    }
    assert(nItemsByPoint.size() == nLocalPoints);

    // ...then compute the global starting indices.
    // Note that we only get the global starting indices for our 
    // own gridpoints.
    auto globalStartingIndices = findGlobalStartingIndices(comm,
                                    commRank,
                                    commSize,
                                    nItemsByPoint);
    assert(globalStartingIndices.size() == nLocalPoints);
    
    // Let everyone know the total number of items.
    // Needed to create the dataset for the actual data.
    // Before the broadcast, globalNumItems is correct only 
    // in the last MPI rank.
    uint32_t globalNumItems = globalStartingIndices.back() + 
                                nItemsByPoint.back();
    MPI_Bcast(&globalNumItems, 1, MPI_UNSIGNED, (commSize - 1), comm);


    // Write the global starting index data to our group.
    // There are nGridPoints values, plus 1 so that we can
    // use the ((i+1) - i) approach to finding how many items per
    // grid point.
    //
    // First, create a 1D dataspace of the right size for the global data...
    SimpleDataSpace<1>::Dimensions globalDims { globalNumItems + 1 };
    SimpleDataSpace<1> startIndicesDataSpace(globalDims);

    // ...then create a dataset for the index data...
    DataSet<FlatStartingIndicesType> startIndicesDataset(*this,
                                        startIndicesDatasetName,
                                        startIndicesDataSpace);

    // ...then specify the shape of the starting index data within
    // our memory space...
    SimpleDataSpace<1>::Dimensions myDims { nLocalPoints };
    SimpleDataSpace<1> memSpace(myDims);

    SimpleDataSpace<1>::Dimensions myOffset { baseX };
    if(commRank == (commSize - 1)) {
        // We're the last rank.  Add an extra element that 
        // contains the total number of data items so that
        // we can use the (index[i+1] - index[i]) approach to
        // determining the number of items per grid point.
        ++myOffset[0];
        globalStartingIndices.push_back(globalNumItems);
    }

    // Select our hyperslab within the file.
    SimpleDataSpace<1> fileSpace(startIndicesDataset);
    H5Sselect_hyperslab(fileSpace.getId(),
                        H5S_SELECT_SET,
                        myOffset.data(),
                        nullptr,
                        myDims.data(),
                        nullptr);

    // Write the data using a collective write.
    PropertyList<H5P_DATASET_XFER> plist;
    H5Pset_dxpl_mpio(plist.getId(), H5FD_MPIO_COLLECTIVE);
    H5Dwrite(startIndicesDatset.getId(),
                TypeInMemory<globalStartIndices::value_type>,
                memSpace.getId(),
                fileSpace.getId(),
                plist.getId(),
                globalStartingIndices.data());
}

HDF5File::Group::FlatStartingIndicesType
HDF5File::Group::findGlobalStartingIndices(MPI_Comm comm,
                        int commRank,
                        int commSize,
                        const FlatStartingIndicesType& localNumItems) const {

    // Second, figure out starting indices for our data items.
    // Note: we're using C++11's std::partial_sum which does an inclusive
    // scan, so we need to manually shift the values in the output
    // vector, and then drop the last value after the scan.
    // If we had C++17 support, we could just use std::exclusive_scan.
    FlatStartingIndicesType localStartingIndices(nLocalGridpoints+1);
    localStartingIndices[0] = 0;
    std::partial_sum(nLocalConcsByGridpoint.begin(),
                        nLocalConcsByGridpoint.end(),
                        localStartingIndices.begin() + 1);
    localStartingIndices.resize(nLocalGridpoints);

    // Third, adjust the local starting indices to make them global
    // starting indices.
    // We do this by 
    uint32_t nLocalConcs = localStartingIndices.back() + 
                                nLocalConcsByGridpoint.back();
    uint32_t globalBaseIdx = 0;
    MPI_Exscan(&nLocalConcs,
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

    FlatStartingIndicesType globalStartingIndices(localStartingIndices.size());
    std::transform(localStartingIndices.begin(), localStartingIndices.end(),
                    globalStartingIndices.begin(),
                    [globalBaseIdx](uint32_t idx) -> uint32_t {
                        return idx + globalBaseIdx;
                    });

    return globalStartingIndices;
}
#endif // READY

} // namespace xolotlCore
#endif // XCORE_HDF5FILE_GROUP_H
