#pragma once

#include <string>

namespace xolotl
{
namespace test
{
class SystemTestCase
{
	static constexpr double defaultTolerance = 1.0e-10;

public:
	SystemTestCase() = delete;

	SystemTestCase(const std::string& caseName,
		const std::string& outputFileName);

	explicit SystemTestCase(const std::string& caseName);

    SystemTestCase&
    tolerance(double tol)
    {
        _tolerance = tol;
        return *this;
    }

	void
	run() const;

	static void
	copyFile(const std::string& fileName);

private:
	bool
	runXolotl() const;

	void
	checkOutput(const std::string& outputFileName,
		const std::string& expectedOutputFileName) const;

private:
	static const std::string _dataDir;
	static const std::string _binDir;
	static const std::string _defaultOutputFileName;

	const std::string _caseName;
	const std::string _outputFileName;

	double _tolerance{defaultTolerance};
};
} // namespace test
} // namespace xolotl
