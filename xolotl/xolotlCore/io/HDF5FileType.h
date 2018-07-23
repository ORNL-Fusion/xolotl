#ifndef XCORE_HDF5FILE_TYPE_H
#define XCORE_HDF5FILE_TYPE_H

#include "xolotlCore/io/HDF5Exception.h"

namespace xolotlCore
{


//----------------------------------------------------------------------------
// int32_t
//----------------------------------------------------------------------------

template<>
inline
HDF5File::TypeInFile<int32_t>::TypeInFile(void)
  : HDF5File::TypeBase("int32_t", H5T_STD_I32LE, false)
{ }

template<>
inline
HDF5File::TypeInMemory<int32_t>::TypeInMemory(void)
  : HDF5File::TypeBase("int32_t", H5T_NATIVE_INT, false)
{ }


//----------------------------------------------------------------------------
// uint32_t
//----------------------------------------------------------------------------

template<>
inline
HDF5File::TypeInFile<uint32_t>::TypeInFile(void)
  : HDF5File::TypeBase("uint32_t", H5T_STD_U32LE, false)
{ }

template<>
inline
HDF5File::TypeInMemory<uint32_t>::TypeInMemory(void)
  : HDF5File::TypeBase("uint32_t", H5T_NATIVE_UINT, false)
{ }


//----------------------------------------------------------------------------
// uint64_t
//----------------------------------------------------------------------------

template<>
inline
HDF5File::TypeInFile<uint64_t>::TypeInFile(void)
  : HDF5File::TypeBase("uint64_t", H5T_STD_U64LE, false)
{ }

template<>
inline
HDF5File::TypeInMemory<uint64_t>::TypeInMemory(void)
  : HDF5File::TypeBase("uint64_t", H5T_NATIVE_ULONG, false)
{ }


#if defined(__clang__) && defined(__APPLE__)
//----------------------------------------------------------------------------
// size_t
//----------------------------------------------------------------------------
// TODO why do I need size_t specialization in addition to uin64_t
// specialization?

template<>
inline
HDF5File::TypeInFile<size_t>::TypeInFile(void)
  : HDF5File::TypeBase("size_t", H5T_STD_U64LE, false)
{ }

template<>
inline
HDF5File::TypeInMemory<size_t>::TypeInMemory(void)
  : HDF5File::TypeBase("size_t", H5T_NATIVE_ULONG, false)
{ }
#endif // defined(__clang__) && defined(__APPLE__)


//----------------------------------------------------------------------------
// std::string
//----------------------------------------------------------------------------

template<>
HDF5File::TypeInMemory<std::string>::TypeInMemory(void);

template<>
HDF5File::TypeInFile<std::string>::TypeInFile(void);


//----------------------------------------------------------------------------
// std::vector<T>
//----------------------------------------------------------------------------


template<typename T>
inline
HDF5File::TypeInFile<std::vector<T> >::TypeInFile(void)
  : HDF5File::TypeBase("std::vector<T>",
                            H5Tvlen_create(TypeInFile<T>().GetId()),
                            true)
{
    if(id < 0)
    {
        throw HDF5Exception("Unable to construct variable-length vector data type for template type T");
    }
}

template<typename T>
inline
HDF5File::TypeInMemory<std::vector<T> >::TypeInMemory(void)
  : HDF5File::TypeBase("std::vector<T>",
                            H5Tvlen_create(TypeInMemory<T>().GetId()),
                            true)
{
    if(id < 0)
    {
        throw HDF5Exception("Unable to construct variable-length vector data type for template type T");
    }
}

} // namespace xolotlCore

#endif // XCORE_HDF5FILE_TYPE_H
