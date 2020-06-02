#ifndef XSOLVER_RANDOMNUMBERGENERATOR_H
#define XSOLVER_RANDOMNUMBERGENERATOR_H

#include <array>

namespace xolotl {
namespace solver {
namespace handler {

// Simple class to wrap a random number generator.
// Intended mainly to allow us to retain the seed used and last few
// values, to support reproducibility.
// But also helps encapsulate the RNG used, so that we can switch
// to a different one if the default rand() proves insufficient.
template<typename ValueType, typename SeedType>
class RandomNumberGenerator {
private:
	static constexpr uint32_t nValuesToSave = 3;

public:
	using RecentValuesType = std::array<ValueType, nValuesToSave>;

private:
	SeedType seed;
	RecentValuesType recentValues;

	ValueType GetRandomValue(void) {
		// Shift recent values to make space for new sample.
		// Be sure to shift from the right so we don't just
		// keep overwriting with the value at slot 0.
		for (auto i = (nValuesToSave - 1); i != 0; --i) {
			recentValues[i] = recentValues[i - 1];
		}

		// Get new sample, save it, and return it to caller.
		recentValues[0] = rand();
		return recentValues[0];
	}

	void SeedGenerator(void) const {
		std::srand(seed);
	}

public:
	RandomNumberGenerator(void) = delete;

	RandomNumberGenerator(SeedType _seed) :
			seed(_seed) {
		SeedGenerator();

		for (auto i = 0; i < nValuesToSave; ++i) {
			recentValues[i] = 0;
		}
	}

	double GetRandomDouble(void) {
		return ((double) GetRandomValue()) / RAND_MAX;
	}

	const RecentValuesType& GetRecentValues(void) const {
		return recentValues;
	}
};

} /* namespace handler */
} /* namespace solver */
} /* namespace xolotl */

#endif // XSOLVER_RANDOMNUMBERGENERATOR_H
