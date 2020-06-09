#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <xolotl/util/Filesystem.h>

namespace testUtils
{
const std::string dataDir = TEST_DATA_DIR;
const std::string binDir = BUILD_DIR;

constexpr double defaultTolerance = 1.0e-10;
const std::string defaultFileName = "retentionOut.txt";

double tolerance = defaultTolerance;

bool
runXolotl(const std::string& caseName)
{
    auto exec = binDir + "/xolotl/xolotl";
    auto paramsFileName = dataDir + "/params_" + caseName + ".txt";
    auto consoleFileName = binDir + "/test/system/cout_" + caseName + ".txt";
    auto command = exec + " " + paramsFileName + " > " + consoleFileName;
    int retCode = std::system(command.c_str());
    return (retCode == 0);
}

bool
copyFile(const std::string& paramsFileName)
{
    int retCode = std::system(("cp " + dataDir + "/" + paramsFileName + " .").c_str());
    return (retCode == 0);
}

std::vector<double>
readOutputFile(const std::string& fileName)
{
    std::ifstream ifs(fileName);
    if (!ifs) {
        throw std::runtime_error("Unable to open file: " + fileName);
    }
    std::string line, tmpStr;
    std::vector<double> ret;
    while (getline(ifs, line)) {
        std::stringstream(line) >> tmpStr;
        if (tmpStr[0] == '#') {
            continue;
        }

        std::stringstream ss(line);
        while (ss >> tmpStr) {
            auto tmp = atof(tmpStr.c_str());
            ret.push_back(tmp);
        }
    }

    return ret;
}

double
diff2Norm(const std::vector<double>& data,
    const std::vector<double>& expectedData)
{
    double diffNorm = 0.0;
    double expectNorm = 0.0;
    auto sz = data.size();
    for (std::size_t i = 0; i < sz; ++i) {
        auto diff = data[i] - expectedData[i];
        auto expect = expectedData[i];
        if (std::isnan(expectedData[i])) {
            expect = 0.0;
            if (std::isnan(data[i])) {
                diff = 0.0;
            }
        }
        if (std::isinf(expectedData[i])) {
            expect = 0.0;
            if (std::isinf(data[i])) {
                diff = 0.0;
            }
        }
        diffNorm += diff*diff;
        expectNorm += expect*expect;
    }
    diffNorm = std::sqrt(diffNorm);
    expectNorm = std::sqrt(expectNorm);
    diffNorm /= expectNorm;
    return diffNorm;
}

double
diffInfNorm(const std::vector<double>& data,
    const std::vector<double>& expectedData)
{
    double diffNorm = 0.0;
    double expectNorm = 0.0;
    auto sz = data.size();
    for (std::size_t i = 0; i < sz; ++i) {
        auto diff = std::abs(data[i] - expectedData[i]);
        auto expect = std::abs(expectedData[i]);
        if (std::isnan(expectedData[i])) {
            expect = 0.0;
            if (std::isnan(data[i])) {
                diff = 0.0;
            }
        }
        if (std::isinf(expectedData[i])) {
            expect = 0.0;
            if (std::isinf(data[i])) {
                diff = 0.0;
            }
        }
        diffNorm = std::max(diffNorm, diff);
        expectNorm = std::max(expectNorm, expect);
    }
    diffNorm /= expectNorm;
    return diffNorm;
}

double
computeDiffNorm(const std::vector<double>& data,
    const std::vector<double>& expectedData)
{
    return diff2Norm(data, expectedData);
}

void
checkOutput(const std::string& outputFileName,
    const std::string& expectedOutputFileName)
{
    auto expectedData = readOutputFile(expectedOutputFileName);
    auto data = readOutputFile(outputFileName);
    BOOST_REQUIRE(expectedData.size() == data.size());
    auto diffNorm = computeDiffNorm(data, expectedData);
    std::cout << std::scientific << std::setprecision(12) << diffNorm << " < " << tolerance << std::endl;
    BOOST_REQUIRE(diffNorm < tolerance);
}

void
runSystemTestCase(const std::string& caseName, const std::string& fileName = defaultFileName,
		double tol = defaultTolerance)
{
    tolerance = tol;
    BOOST_REQUIRE(runXolotl(caseName));

    auto argc = boost::unit_test::framework::master_test_suite().argc;
    auto argv = boost::unit_test::framework::master_test_suite().argv;
    if (argc == 2 && std::strcmp(argv[1], "--approve") == 0) {
        fs::copy_file("./" + fileName,
            dataDir + "/output/" + caseName + ".txt",
            fs::copy_option::overwrite_if_exists);
    }
    else {
        checkOutput("./" + fileName,
            dataDir + "/output/" + caseName + ".txt");
    }
}
}
