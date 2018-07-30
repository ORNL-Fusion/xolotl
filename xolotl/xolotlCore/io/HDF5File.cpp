#include <sstream>
#include <array>
#include "hdf5.h"
#include "xolotlCore/io/HDF5File.h"
#include "xolotlCore/io/HDF5Exception.h"


namespace xolotlCore {

unsigned int
HDF5File::toHDF5AccessMode(AccessMode mode) {

    unsigned int ret = 0;

    switch(mode) {
    case AccessMode::OpenReadOnly:
        ret = H5F_ACC_RDONLY;
        break;

    case AccessMode::OpenReadWrite:
        ret = H5F_ACC_RDWR;
        break;

    case AccessMode::CreateOrTruncateIfExists:
        ret = H5F_ACC_TRUNC;
        break;

    case AccessMode::CreateOrFailIfExists:
        ret = H5F_ACC_EXCL;
        break;

    default:
        throw HDF5Exception("Unrecognized file access mode given when opening/creating file");

    }
    return ret;
}


herr_t
HDF5File::BuildCurrStackLevelErrorString(uint32_t n,
                                        const H5E_error2_t* eptr,
                                        void* client_data)
{
    assert(client_data != nullptr);
    auto& estr = *(static_cast<std::ostringstream*>(client_data));

    estr << n << ':';

    auto constexpr maxMsgSize = 64;
    std::array<char, maxMsgSize> eclass;
    auto ok = H5Eget_class_name(eptr->cls_id, eclass.data(), maxMsgSize);
    estr << ((ok >= 0) ? eclass.data() : "Unknown") << ':';

    std::array<char, maxMsgSize> emajor;
    ok = H5Eget_msg(eptr->maj_num, nullptr, emajor.data(), maxMsgSize);
    estr << ((ok >= 0) ? emajor.data() : "Unknown") << ':';

    std::array<char, maxMsgSize> eminor;
    ok = H5Eget_msg(eptr->min_num, nullptr, eminor.data(), maxMsgSize);
    estr << ((ok >= 0) ? eminor.data() : "Unknown") << '\n';

    return 0;
}


std::string
HDF5File::BuildHDF5ErrorString(void) {

    std::ostringstream estr;

    auto ok = H5Ewalk(H5E_DEFAULT,
                        H5E_WALK_UPWARD,
                        BuildCurrStackLevelErrorString,
                        &estr);
    if(ok < 0)
    {
        estr << "HDF5 error occurred, no details available.";
    }
    return estr.str();
}


HDF5File::HDF5File(fs::path path,
                    AccessMode mode,
                    MPI_Comm _comm,
                    bool par)
  : HDF5Object("/"),
    comm(_comm)
{
    // Obtain the HDF5 flag that corresponds to requested access mode.
    auto hdf5Mode = toHDF5AccessMode(mode);

    // If requested, specify that we want to access the file using parallel I/O.
    PropertyList plist(H5P_FILE_ACCESS);
    auto plistId = plist.getId();
    if(par) {
        // Use parallel I/O for accessing this file.
        H5Pset_fapl_mpio(plistId, comm, MPI_INFO_NULL);
    }
    else {
        // Do not use parallel I/O when accessing this file.
        plistId = H5P_DEFAULT;
    }

    if((mode == AccessMode::OpenReadOnly) or 
            (mode == AccessMode::OpenReadWrite)) {

        // We are opening an existing file.
        setId(H5Fopen(path.string().c_str(), hdf5Mode, plistId));
    }
    else {
        // We are creating/truncating a file.
        setId(H5Fcreate(path.string().c_str(),
                        hdf5Mode,
                        H5P_DEFAULT,    // create plist
                        plistId));       // access mode plist
    }

    if(getId() < 0) {
        throw HDF5Exception(BuildHDF5ErrorString());
    }
}

HDF5File::~HDF5File(void) {
    H5Fclose(getId());
}


bool
HDF5File::hasGroup(fs::path path) const {

    auto cret = H5Lexists(getId(), path.string().c_str(), H5P_DEFAULT);
    if(cret < 0) {
        throw HDF5Exception(BuildHDF5ErrorString());
    }
    return (cret != 0);
}

} // namespace xolotlCore

