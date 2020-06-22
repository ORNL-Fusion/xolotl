#ifndef XCORE_XFILE_TYPE_H
#define XCORE_XFILE_TYPE_H

namespace xolotl
{
namespace io
{
//----------------------------------------------------------------------------
// XFile::TimestepGroup::ConcType
//----------------------------------------------------------------------------

template <>
inline hid_t
HDF5File::TypeInMemory<XFile::TimestepGroup::ConcType>::BuildCompoundType(
	void) const
{
	using ConcType = XFile::TimestepGroup::ConcType;

	auto tid = H5Tcreate(H5T_COMPOUND, sizeof(ConcType));
	H5Tinsert(tid, "ConcType.first", HOFFSET(ConcType, first), H5T_NATIVE_INT);
	H5Tinsert(
		tid, "ConcType.second", HOFFSET(ConcType, second), H5T_NATIVE_DOUBLE);

	return tid;
}

template <>
inline hid_t
HDF5File::TypeInFile<XFile::TimestepGroup::ConcType>::BuildCompoundType(
	void) const
{
	using ConcType = XFile::TimestepGroup::ConcType;

	auto tid = H5Tcreate(H5T_COMPOUND, sizeof(ConcType));
	H5Tinsert(tid, "ConcType.first", HOFFSET(ConcType, first), H5T_STD_I32LE);
	H5Tinsert(
		tid, "ConcType.second", HOFFSET(ConcType, second), H5T_IEEE_F64LE);

	return tid;
}

template <>
inline HDF5File::TypeInFile<XFile::TimestepGroup::ConcType>::TypeInFile(void) :
	HDF5File::TypeBase(
		"XFile::TimestepGroup::ConcType", BuildCompoundType(), true)
{
}

template <>
inline HDF5File::TypeInMemory<XFile::TimestepGroup::ConcType>::TypeInMemory(
	void) :
	HDF5File::TypeBase(
		"XFile::TimestepGroup::ConcType", BuildCompoundType(), true)
{
}

} // namespace io
} // namespace xolotl

#endif // XCORE_XFILE_TYPE_H
