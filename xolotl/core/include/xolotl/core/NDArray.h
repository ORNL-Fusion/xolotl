#ifndef XCORE_NDARRAY_H
#define XCORE_NDARRAY_H

#include <iostream>
#include <sstream>
#include <iterator>
#include <array>

namespace xolotl {
namespace core {

// Forward declaration, apparently needed for the rest to
// compile successfully with some compilers (e.g., clang++ 6).
template<typename T, uint32_t N0, uint32_t ... Ns>
class Array;

// An N-dimensional array of Ts.
// Specialization for 2+ dimensions.
template<typename T, uint32_t N0, uint32_t N1, uint32_t ... Ns>
class Array<T, N0, N1, Ns...> : public std::array<Array<T, N1, Ns...>, N0> {
public:
	Array(void) = default;
	Array(const Array<T, N0, N1, Ns...>& other) = default;

	void Init(const T& val) {
		for (auto i = 0; i < N0; ++i) {
			this->operator[](i).Init(val);
		}
	}

	Array<T, N0, N1, Ns...>& operator+=(const Array<T, N0, N1, Ns...>& other) {
		for (auto i = 0; i < N0; ++i) {
			this->operator[](i) += other[i];
		}
		return *this;
	}

	void DumpTo(std::ostream& os, size_t level = 0) {
		std::ostringstream indentStr;
		for (auto j = 0; j < level; ++j) {
			indentStr << ' ';
		}

		os << indentStr.str() << "[\n";
		for (auto i = 0; i < N0; ++i) {
			this->operator[](i).DumpTo(os, level + 1);
		}
		os << indentStr.str() << "]\n";
	}
};

template<typename T, uint32_t N0, uint32_t N1, uint32_t ...Ns>
Array<T, N0, N1, Ns...> operator+(const Array<T, N0, N1, Ns...>& a,
		const Array<T, N0, N1, Ns...>& b) {
	Array<T, N0, N1, Ns...> ret = a;
	ret += b;
	return ret;
}

// An N-dimensional array of type T.
// Specialization for one dimension.
template<typename T, uint32_t N0>
class Array<T, N0> : public std::array<T, N0> {
public:
	Array(void) = default;
	Array(const Array<T, N0>& other) = default;

	void Init(const T& val) {
		std::fill(this->begin(), this->end(), val);
	}

	Array<T, N0>& operator+=(const Array<T, N0>& other) {
		for (auto i = 0; i < N0; ++i) {
			this->operator[](i) += other[i];
		}
		return *this;
	}

	void DumpTo(std::ostream& os, size_t level = 0) {
		for (auto j = 0; j < level; ++j) {
			os << ' ';
		}
		std::copy(this->begin(), this->end(),
				std::ostream_iterator<T>(os, " "));
		os << std::endl;
	}
};

template<typename T, uint32_t N0>
Array<T, N0> operator+(const Array<T, N0>& a, const Array<T, N0>& b) {
	Array<T, N0> ret = a;
	ret += b;
	return ret;
}

} /* end namespace core */
} /* end namespace xolotl */

#endif // XCORE_NDARRAY_H
