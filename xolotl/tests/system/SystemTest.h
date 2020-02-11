#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <Filesystem.h>

namespace testUtils
{
const std::string dataDir = TEST_DATA_DIR;
const std::string binDir = BUILD_DIR;

constexpr double defaultTolerance = 1.0e-10;

double tolerance = defaultTolerance;

bool
runXolotl(const std::string& paramsFileName)
{
    int retCode = std::system((binDir + "/xolotl " + paramsFileName).c_str());
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

void
checkOutput(const std::string& outputFileName,
    const std::string& expectedOutputFileName)
{
    auto expectedData = readOutputFile(expectedOutputFileName);
    auto data = readOutputFile(outputFileName);
    auto sz = expectedData.size();
    BOOST_REQUIRE(sz == data.size());
    double diffNorm = 0.0;
    for (std::size_t i = 0; i < sz; ++i) {
        auto diff = std::abs(data[i] - expectedData[i]);
        if (std::isnan(data[i]) && std::isnan(expectedData[i])) {
            diff = 0.0;
        }
        if (std::isinf(data[i]) && std::isinf(expectedData[i])) {
            diff = 0.0;
        }
        diffNorm = std::max(diffNorm, diff);
    }
    BOOST_REQUIRE(diffNorm < tolerance);
}

void
runSystemTestCase(const std::string& caseName, double tol = defaultTolerance)
{
    tolerance = tol;
    BOOST_REQUIRE(runXolotl(dataDir + "/params_" + caseName + ".txt"));

    auto argc = boost::unit_test::framework::master_test_suite().argc;
    auto argv = boost::unit_test::framework::master_test_suite().argv;
    if (argc == 2 && std::strcmp(argv[1], "--approve") == 0) {
        fs::copy_file("./retentionOut.txt",
            dataDir + "/output/retention_" + caseName + ".txt",
            fs::copy_option::overwrite_if_exists);
    }
    else {
        checkOutput("./retentionOut.txt",
            dataDir + "/output/retention_" + caseName + ".txt");
    }
}
}
