#ifndef XCORE_INTEGER_RANGE_H
#define XCORE_INTEGER_RANGE_H

#include <tuple>
#include <stdexcept>

namespace xolotl {
namespace core {

/**
 * A range of integers [l, h).
 * Only the bounds are explicitly stored.
 * Implementation based on responses to posted question on StackOverflow:
 * https://stackoverflow.com/questions/7185437/is-there-a-range-class-in-c11-for-use-with-range-based-for-loops
 */
template<typename T, typename = typename std::enable_if<
		std::is_integral<T>::value, T>::type>
class IntegerRange: public std::pair<T, T> {
private:
	/**
	 * Our range endpoints.  Range is closed on lower side, open on higher side.
	 */
	T lo;
	T hi;

public:
	/**
	 * An iterator to values within our range.
	 */
	struct iterator {
		T val;

		iterator(const T& _val) :
				val(_val) {
			// Nothing else to do.
		}

		const iterator& operator++() {
			++val;
			return *this;
		}

		bool operator!=(const iterator& _other) {
			return val != _other.val;
		}

		T operator*() const {
			return val;
		}
	};

	/**
	 * Default ctor, disallowed because we must construct with
	 * values for the bounds.
	 */
	IntegerRange() = delete;

	/**
	 * Construct an integer range.
	 */
	IntegerRange(const T& _lo, const T& _hi) :
			lo(_lo), hi(_hi) {
		// Verify that lo is not larger than hi.
		if (lo > hi) {
			throw std::out_of_range("lower bound greater than higher bound");
		}
	}

	/**
	 * Copy another range.
	 */
	IntegerRange(const IntegerRange& _other) :
			lo(_other.lo), hi(_other.hi) {
		// Nothing else to do.
	}

	/**
	 * Obtain iterator to beginning of range.
	 *
	 * @return Iterator to first value in range.
	 */
	iterator begin(void) const {
		return lo;
	}

	/**
	 * Obtain iterator to end of range.
	 * Follows C++ convention of end being one past last valid value.
	 *
	 * @return Iterator to one past last value in range.
	 */
	iterator end(void) const {
		return hi;
	}

	/**
	 * Check whether our range contains the given value.
	 *
	 * @param _val The value to check whether it is in our range.
	 * @return True iff _val is within our range.
	 */
	bool contains(const T& _val) const {
		return ((_val >= lo) and (_val < hi));
	}
};

} /* end namespace core */
} /* end namespace xolotl */

#endif // XCORE_INTEGER_RANGE_H
