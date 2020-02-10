#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

const std::string dataDir = TEST_DATA_DIR;
const std::string binDir = BUILD_DIR;

constexpr double tolerance = 1.0e-10;

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
    std::string tmpStr;
    double tmp;
    std::vector<double> ret;
    while (getline(ifs, tmpStr)) {
        std::string dummy;
        std::stringstream(tmpStr) >> dummy;
        if (dummy[0] == '#') {
            continue;
        }

        std::stringstream ss(tmpStr);
        while (ss >> tmp) {
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
        diffNorm = std::max(diffNorm, diff);
    }
    BOOST_REQUIRE(diffNorm < tolerance);
}

void
runSystemTestCase(const std::string& caseName)
{
    BOOST_REQUIRE(runXolotl(dataDir + "/params_" + caseName + ".txt"));
    checkOutput("./retentionOut.txt",
        dataDir + "/output/retention_" + caseName + ".txt");
}

BOOST_AUTO_TEST_CASE(NE_9999DOF_plsm) {
    runSystemTestCase("NE_9999DOF_plsm");
}
