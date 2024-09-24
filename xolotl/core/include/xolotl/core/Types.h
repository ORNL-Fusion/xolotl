#pragma once

#include <type_traits>

#include <Kokkos_Core.hpp>
#include <Kokkos_OffsetView.hpp>

#include <xolotl/config.h>
#include <xolotl/util/Array.h>

namespace xolotl
{
namespace core
{
using RowColPair = util::Array<IdType, 2>;

template <typename T>
class SizedArrayWrapper
{
public:
	// using reference = T&;
	// using const_reference = std::add_const_t<T>&;
	// using size_type = size_t;
	// using difference_type = ptrdiff_t;
	// using value_type = T;
	// using pointer = T*;
	// using const_pointer = std::add_const_t<T>*;

	KOKKOS_INLINE_FUNCTION
	SizedArrayWrapper(T* data, std::size_t size) : _data(data), _size(size)
	{
	}

	KOKKOS_INLINE_FUNCTION
	SizedArrayWrapper&
	operator=(const SizedArrayWrapper<T>& rhs)
	{
		const size_t n = size() < rhs.size() ? size() : rhs.size();
		for (size_t i = 0; i < n; ++i) {
			_data[i] = rhs[i];
		}
		return *this;
	}

	template <typename TIndex>
	KOKKOS_INLINE_FUNCTION
	T&
	operator[](const TIndex& i)
	{
		static_assert((std::is_integral_v<TIndex> || std::is_enum_v<TIndex>),
			"Must be integral argument");
		assert(i < _size);
		return _data[i];
	}

	template <typename TIndex>
	KOKKOS_INLINE_FUNCTION
	const T&
	operator[](const TIndex& i) const
	{
		static_assert((std::is_integral_v<TIndex> || std::is_enum_v<TIndex>),
			"Must be integral argument");
		assert(i < _size);
		return _data[i];
	}

	KOKKOS_INLINE_FUNCTION
	std::size_t
	size() const
	{
		return _size;
	}

	KOKKOS_INLINE_FUNCTION
	bool
	empty() const
	{
		return 0 == _size;
	}

	KOKKOS_INLINE_FUNCTION
	std::size_t
	max_size() const
	{
		return _size;
	}

	KOKKOS_INLINE_FUNCTION
	T*
	data()
	{
		return _data;
	}

	KOKKOS_INLINE_FUNCTION
	const T*
	data() const
	{
		return _data;
	}

private:
	T* _data;
	std::size_t _size;
};

using StencilConcArray = SizedArrayWrapper<Kokkos::View<const double*>>;
} // namespace core

using DefaultMemSpace = Kokkos::DefaultExecutionSpace::memory_space;
template <typename T>
using PetscOffsetView =
	Kokkos::Experimental::OffsetView<T, Kokkos::LayoutRight, DefaultMemSpace>;

} // namespace xolotl
