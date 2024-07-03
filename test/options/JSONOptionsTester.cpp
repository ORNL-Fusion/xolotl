#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/framework.hpp>
#include <boost/test/unit_test.hpp>

#include <xolotl/options/InvalidOptionValue.h>
#include <xolotl/options/JSONOptions.h>
#include <xolotl/options/Options.h>
#include <xolotl/test/config.h>
const std::string dataDir = TO_STRING(XOLOTL_TEST_DATA_DIR);

using namespace std;
using namespace xolotl;
using namespace xolotl::options;

/**
 * Test suite for the Options class.
 */
BOOST_AUTO_TEST_SUITE(JSONOptions_testSuite)

BOOST_AUTO_TEST_CASE(noOptions)
{
	JSONOptions opts;

	// Build a fake, empty command line.
	// Here, "empty" is what the shell would give us if no arguments
	// were provided (i.e., it contains the executable name or path),
	// and we skipped that executable name before calling the
	// parsing method as specified in its comment.
	const char* argv[] = {"./xolotl"};

	// Attempt to read the parameters
	BOOST_CHECK_THROW(opts.readParams(1, argv), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(badParamFileName)
{
	JSONOptions opts;

	std::string fileName("bla.json");

	// Build a command line with a non existing parameter file
	const char* argv[] = {"./xolotl", fileName.c_str()};

	// Attempt to read the parameters
	BOOST_CHECK_THROW(opts.readParams(2, argv), std::runtime_error);
}

// TODO: currently don't do anything for unexpected parameters
// BOOST_AUTO_TEST_CASE(badParamFile)
// {
// 	try {
// 		JSONOptions opts;

// 		// Create a bad parameter file
//         std::string fileName = "param_bad.json";
// 		std::ofstream badParamFile(fileName);
// 		badParamFile << "{ \"netFile\": \"tungsten.txt\" }\n";
// 		badParamFile.close();

// 		// Build a command line with a parameter file containing bad options
// 		const char* argv[] = {"./xolotl", fileName.c_str()};

// 		// Attempt to read the parameter file
// 		opts.readParams(2, argv);

// 		// Remove the created file
// 		fs::remove(fileName);
// 	}
// 	catch (const std::exception& e) {
// 		// Great
// 		std::cerr << e.what() << std::endl;
// 		return;
// 	}

// 	// Should not get here
// 	BOOST_REQUIRE(false);
// }

BOOST_AUTO_TEST_CASE(goodParamFile)
{
	JSONOptions opts;

	// Create a good parameter file
	std::string fileName = "param_good.json";
	std::ofstream goodParamFile(fileName);
	goodParamFile
		<< "{\n"
		<< "\"vizHandler\": \"std\",\n"
		<< "\"petscArgs\": [\n"
		<< "\"-fieldsplit_0_pc_type redundant\",\n"
		<< "\"-ts_max_snes_failures 200\",\n"
		<< "\"-pc_fieldsplit_detect_coupling\",\n"
		<< "\"-ts_adapt_dt_max 10\",\n"
		<< "\"-pc_type fieldsplit\",\n"
		<< "\"-fieldsplit_1_pc_type sor\",\n"
		<< "\"-ts_final_time 1000\",\n"
		<< "\"-ts_max_steps 3\"\n"
		<< "],\n"
		<< "\"checkpointFile\": \"tungsten.txt\",\n"
		<< "\"tempHandler\": \"constant\",\n"
		<< "\"tempParam\": 900,\n"
		<< "\"perfHandler\": \"os\",\n"
		<< "\"flux\": 1.5,\n"
		<< "\"material\": \"W100\",\n"
		<< "\"initialConc\": \"V 1 0.05\",\n"
		<< "\"dimensions\": 1,\n"
		<< "\"gridType\": \"nonuniform\",\n"
		<< "\"gridParam\": 10,\n"
		<< "\"process\": \"diff\",\n"
		<< "\"grouping\": \"11 2 4\",\n"
		<< "\"sputtering\": 0.5,\n"
		<< "\"boundary\": \"1 1\",\n"
		<< "\"heatLossPortion\": 0.5,\n"
		<< "\"burstingDepth\": 5.0,\n"
		<< "\"burstingFactor\": 2.5,\n"
		<< "\"zeta\": 0.6,\n"
		<< "\"radiusSize\": \"5 0 3\",\n"
		<< "\"density\": 9.0,\n"
		<< "\"lattice\": 0.1,\n"
		<< "\"impurityRadius\": 0.5,\n"
		<< "\"biasFactor\": 2.0,\n"
		<< "\"hydrogenFactor\": 0.5,\n"
		<< "\"xenonDiffusivity\": 3.0,\n"
		<< "\"fissionYield\": 0.3,\n"
		<< "\"heVRatio\": 5.0,\n"
		<< "\"migrationThreshold\": 1.0,\n"
		<< "\"fluxDepthProfileFilePath\": \"path/to/flux/profile/file.txt\",\n"
		<< "\"basalPortion\": 0.6,\n"
		<< "\"transitionSize\": 300,\n"
		<< "\"cascadeDose\": 5000.0,\n"
		<< "\"cascadeEfficiency\": 0.2\n"
		<< "}\n";
	goodParamFile.close();

	// Build a command line with a parameter file containing good options
	const char* argv[] = {"./xolotl", fileName.c_str()};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check the checkpoint filename
	BOOST_REQUIRE_EQUAL(opts.getCheckpointFilePath(), "tungsten.txt");

	// Check the temperature
	BOOST_REQUIRE_EQUAL(opts.getTempHandlerName(), "constant");
	BOOST_REQUIRE_EQUAL(opts.getTempParam(), 900.0);
	BOOST_REQUIRE_EQUAL(opts.getTempParam(1), 0.0);

	// Check if the flux option is used
	BOOST_REQUIRE_EQUAL(opts.useFluxAmplitude(), true);
	BOOST_REQUIRE_EQUAL(opts.getFluxAmplitude(), 1.5);

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getPerfHandlerName(), "os");

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getVizHandlerName(), "std");

	// Check the material option
	BOOST_REQUIRE_EQUAL(opts.getMaterial(), "W100");

	// Check the initial vacancy concentration option
	BOOST_REQUIRE_EQUAL(opts.getInitialConcentration(), "V 1 0.05");

	// Check the number of dimensions option
	BOOST_REQUIRE_EQUAL(opts.getDimensionNumber(), 1);

	// Check the grid options
	BOOST_REQUIRE_EQUAL(opts.getGridTypeName(), "nonuniform");
	BOOST_REQUIRE_EQUAL(opts.getGridParam(0), 10);

	// Check the grouping option
	BOOST_REQUIRE_EQUAL(opts.getGroupingMin(), 11);
	BOOST_REQUIRE_EQUAL(opts.getGroupingWidthA(), 2);
	BOOST_REQUIRE_EQUAL(opts.getGroupingWidthB(), 4);

	// Check the sputtering option
	BOOST_REQUIRE_EQUAL(opts.getSputteringYield(), 0.5);

	// Check the bursting depth option
	BOOST_REQUIRE_EQUAL(opts.getBurstingDepth(), 5.0);

	// Check the bursting factor option
	BOOST_REQUIRE_EQUAL(opts.getBurstingFactor(), 2.5);

	// Check the boundary conditions
	BOOST_REQUIRE_EQUAL(opts.getLeftBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getRightBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getBottomBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getTopBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getFrontBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getBackBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getBCString(), "mirror");
	BOOST_REQUIRE_EQUAL(opts.getHeatLossPortion(), 0.5);

	// Check the electronic stopping power option
	BOOST_REQUIRE_EQUAL(opts.getZeta(), 0.6);

	// Check the minimum size option for the average radius computation
	auto sizes = opts.getRadiusMinSizes();
	BOOST_REQUIRE_EQUAL(sizes[0], 5);
	BOOST_REQUIRE_EQUAL(sizes[1], 0);
	BOOST_REQUIRE_EQUAL(sizes[2], 3);

	// Check the xenon density option
	BOOST_REQUIRE_EQUAL(opts.getDensity(), 9.0);

	// Check the lattice parameter option
	BOOST_REQUIRE_EQUAL(opts.getLatticeParameter(), 0.1);

	// Check the impurity radius option
	BOOST_REQUIRE_EQUAL(opts.getImpurityRadius(), 0.5);

	// Check the interstitial bias factor option
	BOOST_REQUIRE_EQUAL(opts.getBiasFactor(), 2.0);

	// Check the hydrogen radius factor option
	BOOST_REQUIRE_EQUAL(opts.getHydrogenFactor(), 0.5);

	// Check the xenon diffusivity option
	BOOST_REQUIRE_EQUAL(opts.getXenonDiffusivity(), 3.0);

	// Check the fission yield option
	BOOST_REQUIRE_EQUAL(opts.getFissionYield(), 0.3);

	// Check the HeV ration option
	BOOST_REQUIRE_EQUAL(opts.getHeVRatio(), 5.0);

	// Check the migration threshold option
	BOOST_REQUIRE_EQUAL(opts.getMigrationThreshold(), 1.0);

	// Check the flux depth profile filename
	BOOST_REQUIRE_EQUAL(
		opts.getFluxDepthProfileFilePath(), "path/to/flux/profile/file.txt");

	// Check the basal portion
	BOOST_REQUIRE_EQUAL(opts.getBasalPortion(), 0.6);

	// Check the transition size
	BOOST_REQUIRE_EQUAL(opts.getTransitionSize(), 300);

	// Check the cascade options
	BOOST_REQUIRE_EQUAL(opts.getCascadeDose(), 5000);
	BOOST_REQUIRE_EQUAL(opts.getCascadeEfficiency(), 0.2);

	// Check the physical processes option
	auto map = opts.getProcesses();
	BOOST_REQUIRE_EQUAL(map["diff"], true);
	BOOST_REQUIRE_EQUAL(map["advec"], false);
	BOOST_REQUIRE_EQUAL(map["modifiedTM"], false);
	BOOST_REQUIRE_EQUAL(map["movingSurface"], false);
	BOOST_REQUIRE_EQUAL(map["reaction"], false);
	BOOST_REQUIRE_EQUAL(map["attenuation"], false);
	BOOST_REQUIRE_EQUAL(map["bursting"], false);
	BOOST_REQUIRE_EQUAL(map["resolution"], false);
	BOOST_REQUIRE_EQUAL(map["heterogeneous"], false);
	BOOST_REQUIRE_EQUAL(map["constant"], false);

	// Check the PETSc options
	BOOST_REQUIRE_EQUAL(opts.getPetscArg(),
		"-fieldsplit_0_pc_type redundant "
		"-ts_max_snes_failures 200 "
		"-pc_fieldsplit_detect_coupling "
		"-ts_adapt_dt_max 10 -pc_type fieldsplit "
		"-fieldsplit_1_pc_type sor -ts_final_time 1000 "
		"-ts_max_steps 3");

	// Remove the created file
	fs::remove(fileName);
}

BOOST_AUTO_TEST_CASE(goodParamFileNoHDF5)
{
	JSONOptions opts;

	// Create a good parameter file
	std::string fileName = "param_good.json";
	std::ofstream goodParamFile(fileName);
	goodParamFile << "{ \"netParam\": \"8 1 0 5 3\" }\n";
	goodParamFile.close();

	// Build a command line with a parameter file containing good options
	const char* argv[] = {"./xolotl", fileName.c_str()};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check if we use the HDF5 file
	BOOST_REQUIRE_EQUAL(opts.useHDF5(), false);

	// Check the network parameters
	BOOST_REQUIRE_EQUAL(opts.getMaxImpurity(), 8);
	BOOST_REQUIRE_EQUAL(opts.getMaxD(), 1);
	BOOST_REQUIRE_EQUAL(opts.getMaxT(), 0);
	BOOST_REQUIRE_EQUAL(opts.getMaxV(), 5);
	BOOST_REQUIRE_EQUAL(opts.getMaxI(), 3);
	BOOST_REQUIRE_EQUAL(opts.getMaxPureV(), 5);

	// Remove the created file
	fs::remove(fileName);
}

BOOST_AUTO_TEST_CASE(wrongPerfHandler)
{
	JSONOptions opts;

	// Create a parameter file with a wrong performance handler name
	std::string fileName = "param_perf_wrong.json";
	std::ofstream paramFile(fileName);
	paramFile << "{ \"perfHandler\": \"bogus\" }\n";
	paramFile.close();

	// Build a command line with a parameter file containing a wrong performance
	// handler option
	const char* argv[] = {"./xolotl", fileName.c_str()};

	// Attempt to read the parameter file
	BOOST_CHECK_THROW(opts.readParams(2, argv), InvalidOptionValue);

	// Remove the created file
	fs::remove(fileName);
}

BOOST_AUTO_TEST_CASE(wrongVizHandler)
{
	JSONOptions opts;

	// Create a parameter file with a wrong visualization handler name
	std::string fileName = "param_viz_wrong.json";
	std::ofstream paramFile(fileName);
	paramFile << "{ \"vizHandler\": \"bogus\" }\n";
	paramFile.close();

	// Build a command line with a parameter file containing a wrong performance
	// handler option
	const char* argv[] = {"./xolotl", fileName.c_str()};

	// Attempt to read the parameter file
	BOOST_CHECK_THROW(opts.readParams(2, argv), InvalidOptionValue);

	// Remove the created file
	fs::remove(fileName);
}

BOOST_AUTO_TEST_CASE(goodParamFileWithProfiles)
{
	// Create a file with temperature profile data
	// First column with the time and the second with
	// the temperature at that time.
	std::string tempFileName = "temperatureFile.dat";
	std::ofstream writeTempFile(tempFileName);
	writeTempFile << "0.0 2.0 \n"
					 "1.0 1.99219766723 \n"
					 "2.0 1.87758256189 \n"
					 "3.0 1.4311765168 \n"
					 "4.0 0.583853163453 \n"
					 "5.0 0.000137654918313 \n"
					 "6.0 0.789204200569 \n"
					 "7.0 1.9875147713 \n"
					 "8.0 0.854499966191 \n"
					 "9.0 0.235300873168 \n"
					 "10.0 1.99779827918";
	writeTempFile.close();

	// Create a file with a time profile for the flux
	// First column with the time and the second with
	// the amplitude (in He/nm2/s) at that time.
	std::string fluxFileName = "fluxFile.dat";
	std::ofstream writeFluxFile(fluxFileName);
	writeFluxFile << "0.0 1000.0 \n"
					 "1.0 4000.0 \n"
					 "2.0 2000.0 \n"
					 "3.0 3000.0 \n"
					 "4.0 0.0";
	writeFluxFile.close();

	JSONOptions opts;

	// Create a parameter file using these two profile files
	std::string fileName = "param_good_profiles.json";
	std::ofstream paramFile(fileName);
	paramFile << "{\n"
			  << "\"fluxFile\": " << std::quoted(fluxFileName) << ",\n"
			  << "\"tempHandler\": \"profile\",\n"
			  << "\"tempFile\": " << std::quoted(tempFileName) << "\n"
			  << "}\n";
	paramFile.close();

	// Build a command line with a parameter file containing
	// the two profile options
	const char* argv[] = {"./xolotl", fileName.c_str()};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check the temperature
	BOOST_REQUIRE_EQUAL(opts.getTempHandlerName(), "profile");
	BOOST_REQUIRE_EQUAL(opts.getTempProfileFilename(), tempFileName);

	// Check if the heFlux option is used
	BOOST_REQUIRE_EQUAL(opts.useFluxAmplitude(), false);

	// Check if the time profile option is used for the flux
	BOOST_REQUIRE_EQUAL(opts.useFluxTimeProfile(), true);
	BOOST_REQUIRE_EQUAL(opts.getFluxTimeProfileFilePath(), fluxFileName);

	// Remove the created files
	fs::remove(tempFileName);
	fs::remove(fluxFileName);
	fs::remove(fileName);
}

BOOST_AUTO_TEST_CASE(readGridIn)
{
	// Create a grid file
	std::string gridFileName("grid.dat");
	std::ofstream gridFile(gridFileName);
	gridFile << "0.0 1.0 3.0 7.0";
	gridFile.close();

	JSONOptions opts;

	// Create a parameter file with a file name for the grid
	std::string fileName = "param_read_in_grid.json";
	std::ofstream paramFile(fileName);
	paramFile << "{\n"
			  << "\"gridType\": \"read\",\n"
			  << "\"gridFile\": \"grid.dat\"\n"
			  << "}\n";
	paramFile.close();

	// Build a command line with a parameter file containing a wrong performance
	// handler option
	const char* argv[] = {"./xolotl", fileName.c_str()};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check the grid options
	BOOST_REQUIRE_EQUAL(opts.getGridTypeName(), "read");
	BOOST_REQUIRE_EQUAL(opts.getGridFilename(), "grid.dat");

	// Remove the created file
	fs::remove(fileName);
	fs::remove(gridFileName);
}

BOOST_AUTO_TEST_CASE(wrongFluxProfile)
{
	JSONOptions opts;

	// Create a parameter file with a wrong flux profile file name name
	std::string fileName = "param_flux_wrong.json";
	std::ofstream paramFile(fileName);
	paramFile << "{ \"fluxFile\": \"bogus\" }\n";
	paramFile.close();

	// Build a command line with a parameter file containing a wrong performance
	// handler option
	const char* argv[] = {"./xolotl", fileName.c_str()};

	// Attempt to read the parameter file
	BOOST_CHECK_THROW(opts.readParams(2, argv), std::runtime_error);

	// Remove the created file
	fs::remove(fileName);
}

BOOST_AUTO_TEST_CASE(wrongTempProfile)
{
	JSONOptions opts;

	// Create a parameter file with a wrong temperature profile file name name
	std::string fileName = "param_temp_wrong.json";
	std::ofstream paramFile(fileName);
	paramFile << "{ \"tempFile\": \"bogus\" }\n";
	paramFile.close();

	// Build a command line with a parameter file containing a wrong performance
	// handler option
	const char* argv[] = {"./xolotl", fileName.c_str()};

	// Attempt to read the parameter file
	BOOST_CHECK_THROW(opts.readParams(2, argv), std::runtime_error);

	// Remove the created file
	fs::remove(fileName);
}

BOOST_AUTO_TEST_CASE(papiPerfHandler)
{
	JSONOptions opts;

	// Create a parameter file using the PAPI performance handlers
	std::string fileName = "param_good_perf_papi.json";
	std::ofstream paramFile(fileName);
	paramFile << "{ \"perfHandler\": \"papi\" }\n";
	paramFile.close();

	// Build a command line with a parameter file
	const char* argv[] = {"./xolotl", fileName.c_str()};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getPerfHandlerName(), "papi");

	// Remove the created file
	fs::remove(fileName);
}

BOOST_AUTO_TEST_CASE(osPerfHandler)
{
	JSONOptions opts;

	// Create a parameter file using the OS performance handlers
	std::string fileName = "param_good_perf_os.json";
	std::ofstream paramFile(fileName);
	paramFile << "{ \"perfHandler\": \"os\" }\n";
	paramFile.close();

	// Build a command line with a parameter file
	const char* argv[] = {"./xolotl", fileName.c_str()};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getPerfHandlerName(), "os");

	// Remove the created file
	fs::remove(fileName);
}

BOOST_AUTO_TEST_CASE(dummyPerfHandler)
{
	JSONOptions opts;

	// Create a parameter file using the dummy performance handlers
	std::string fileName = "param_good_perf_dummy.json";
	std::ofstream paramFile(fileName);
	paramFile << "{ \"perfHandler\": \"dummy\" }\n";
	paramFile.close();

	// Build a command line with a parameter file containing good options
	const char* argv[] = {"./xolotl", fileName.c_str()};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getPerfHandlerName(), "dummy");

	// Remove the created file
	fs::remove(fileName);
}

BOOST_AUTO_TEST_SUITE_END()
