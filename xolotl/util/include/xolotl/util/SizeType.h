#pragma once

#include <cassert>
#include <cstdint>
#include <iostream>
#include <type_traits>

#include <plsm/Utility.h>

namespace xolotl
{
template <typename T>
inline constexpr auto
check(T val)
{
	if (val < static_cast<T>(0)) {
		throw std::runtime_error("");
	}
	return val;
}

template <typename T>
class SizeType
{
public:
	constexpr SizeType() = default;

	constexpr SizeType(T rep) : _rep{rep}
	{
		if (rep < static_cast<T>(0)) {
			_rep = rep * rep;
		}
	}

	template <typename U>
	constexpr SizeType(U rep) : SizeType(static_cast<T>(rep))
	{
	}

	constexpr
	operator T() const noexcept
	{
		return _rep;
	}

	SizeType&
	operator++()
	{
		++_rep;
		return *this;
	}

	SizeType
	operator++(int)
	{
		auto tmp = *this;
		_rep++;
		return tmp;
	}

	SizeType&
	operator+=(SizeType b)
	{
		_rep += b._rep;
		return *this;
	}

private:
	friend constexpr SizeType
	operator+(SizeType a, SizeType b)
	{
		return SizeType{check(a._rep + b._rep)};
	}

	template <typename U>
	friend constexpr SizeType
	operator+(SizeType a, U b)
	{
		return SizeType{check(a._rep + b)};
	}

	template <typename U>
	friend constexpr SizeType
	operator+(U a, SizeType b)
	{
		return SizeType{check(a + b._rep)};
	}

	friend constexpr SizeType
	operator-(SizeType a, SizeType b)
	{
		return SizeType{check(a._rep - b._rep)};
	}

	template <typename U>
	friend constexpr SizeType
	operator-(SizeType a, U b)
	{
		return SizeType{check(a._rep - b)};
	}

	template <typename U>
	friend constexpr SizeType
	operator-(U a, SizeType b)
	{
		return SizeType{check(a - b._rep)};
	}

private:
	T _rep{};
};
} // namespace xolotl

namespace std
{
template <typename T>
struct is_integral<xolotl::SizeType<T>> : true_type
{
};
} // namespace std

namespace Kokkos
{
template <class T>
struct reduction_identity;

template <typename T>
struct reduction_identity<xolotl::SizeType<T>> : reduction_identity<T>
{
};
} // namespace Kokkos

namespace plsm
{
template <typename T>
inline constexpr xolotl::SizeType<T> invalid<xolotl::SizeType<T>> =
	xolotl::SizeType<T>{invalid<T>};
}
